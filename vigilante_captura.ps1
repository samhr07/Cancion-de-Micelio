# Vigilante de captura_estacional.  Se registra como tarea programada cada 5 min.
#
# El 2026-08-11 la captura murio a las 20:45 y nadie lo vio hasta las 06:41 del
# dia siguiente: 10 h de mercado perdidas por no mirar.
#
# DOS DEFECTOS ENCONTRADOS AL PROBARLO, los dos anotados porque volverian a
# morder a quien reescriba esto:
#
# 1. LA RUTA NO PUEDE IR COMO LITERAL. El directorio del proyecto lleva una "o"
#    acentuada. PowerShell 5.1 lee un .ps1 sin BOM como ANSI, asi que el literal
#    llega mutilado y `Start-Process -WorkingDirectory` revienta con
#    DirectoryNotFoundException. Se usa `$PSScriptRoot`, que lo entrega el host
#    ya decodificado. Todo este archivo es ASCII a proposito.
#
# 2. NO SE DECIDE POR LA TABLA DE PROCESOS. `Get-CimInstance Win32_Process`
#    falla en esta maquina de forma intermitente (error de resolucion de SID), y
#    con `-EA SilentlyContinue` devuelve VACIO -- indistinguible de "no corre".
#    En la primera version eso llevaba directo a arrancar una SEGUNDA instancia
#    sobre el mismo directorio, que se pisarian: la numeracion de partes es
#    `len([f for f in os.listdir(sub) if f.endswith('.parquet')])`, o sea que
#    las dos calcularian el mismo indice. Es la misma familia del defecto de las
#    dos instancias de Micelio.py de la v1.3.
#
#    Se decide por EL DATO, que es lo que de verdad importa y ademas no puede
#    mentir: si hay parquet fresco, hay captura sana, corra quien corra.

$dir     = $PSScriptRoot
$py      = "C:\Users\Usuario\miniconda3\python.exe"
$destino = Join-Path $dir "telemetria\estacional"
$reg     = Join-Path $dir "telemetria\vigilante.log"
$marca   = Join-Path $dir "telemetria\vigilante_ultimo_arranque.txt"

$MIN_SIN_VOLCAR   = 25    # el volcado es cada 10 min; 25 es cuelgue real
$MIN_TRAS_ARRANCAR = 20   # margen para que un recien arrancado llene su buffer

function Registro($msg) {
    "{0}  {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $msg |
        Out-File -FilePath $reg -Append -Encoding utf8
}

# --- 0. Respaldo incremental a DOS destinos, MEJOR ESFUERZO -----------------
# Va primero y entero dentro de try/catch: si un respaldo falla, el vigilante
# tiene que seguir haciendo su trabajo. Un respaldo caido es una molestia; una
# captura caida son horas de mercado.
#
#   D:\Micelio                              proyecto ENTERO, incluido .git.
#                                           Otro medio fisico: es el unico que
#                                           sobrevive a que muera el SSD.
#   C:\Users\Usuario\respaldo_micelio       solo telemetria, mismo disco. Cubre
#                                           el hueco del anterior: que la USB
#                                           este desconectada o se pierda.
#
# `/XO` copia solo lo que no esta ya, y los parquet no se reescriben nunca, asi
# que en regimen son un pu?ado de archivos por pasada. Estrangulado a una vez
# por hora y por destino, para no releer 2500 nombres cada 5 min.
#
# [!] `/FFT` ES OBLIGATORIO HACIA LA USB. Esta formateada en FAT32, que marca
#     tiempos con resolucion de 2 s, contra los 100 ns de NTFS. Sin `/FFT` cada
#     archivo parece distinto por esa diferencia espuria y robocopy vuelve a
#     copiar 1 GB ENTERO cada hora, sobre una USB que da 226 MB/min.
$destinos = @(
    @{ nombre = "USB";  origen = $dir;                          destino = "D:\Micelio";
       marca = "ultimo_respaldo_usb.txt";  fat = $true;  raiz = "D:\" },
    @{ nombre = "disco"; origen = (Join-Path $dir "telemetria"); destino = "C:\Users\Usuario\respaldo_micelio\telemetria";
       marca = "ultimo_respaldo_disco.txt"; fat = $false; raiz = "C:\" }
)

foreach ($dd in $destinos) {
    try {
        $marcaResp = Join-Path $dir ("telemetria\" + $dd.marca)
        $toca = $true
        if (Test-Path $marcaResp) {
            $tr = [datetime]::MinValue
            if ([datetime]::TryParse((Get-Content $marcaResp -Raw).Trim(), [ref]$tr)) {
                $toca = ((Get-Date) - $tr).TotalMinutes -ge 60
            }
        }
        if (-not $toca) { continue }

        # La USB puede no estar puesta. Se comprueba ANTES de invocar robocopy:
        # si no esta, se anota y se sigue -- no es un fallo del vigilante.
        #
        # [!] SE ANOTA EL CAMBIO DE ESTADO, NO CADA PASADA. La primera version
        #     escribia una linea por hora, y la USB estuvo 2.7 dias fuera: 63
        #     lineas identicas que nadie leyo. Un log que repite lo mismo deja
        #     de ser un log. Ahora solo entra la TRANSICION (se fue / volvio),
        #     con lo que el archivo pasa a ser un registro de sucesos.
        $estadoPrev = Join-Path $dir ("telemetria\_disp_" + $dd.nombre + ".txt")
        $hay = Test-Path $dd.raiz
        $antes = if (Test-Path $estadoPrev) { (Get-Content $estadoPrev -Raw).Trim() } else { "" }
        if ("$hay" -ne $antes) {
            Registro ("RESPALDO {0}: {1} {2}" -f $dd.nombre, $dd.raiz,
                      $(if ($hay) { "REAPARECIO. Se sincroniza." } else { "DESAPARECIO. Respaldo detenido hasta que vuelva." }))
            "$hay" | Out-File -FilePath $estadoPrev -Encoding ascii
        }
        if (-not $hay) {
            (Get-Date).ToString("o") | Out-File -FilePath $marcaResp -Encoding ascii
            continue
        }

        $extra = @("/E", "/XO", "/R:1", "/W:2", "/XD", "__pycache__",
                   "/NFL", "/NDL", "/NP", "/NJH", "/NJS", "/MT:4")
        if ($dd.fat) { $extra += "/FFT" }
        $null = robocopy $dd.origen $dd.destino @extra
        # robocopy devuelve 0-7 en exito (1 = se copio algo). >= 8 es error.
        if ($LASTEXITCODE -ge 8) {
            Registro ("RESPALDO {0}: robocopy devolvio {1}" -f $dd.nombre, $LASTEXITCODE)
        } else {
            (Get-Date).ToString("o") | Out-File -FilePath $marcaResp -Encoding ascii
        }
    } catch {
        Registro ("RESPALDO {0} fallido (no bloquea): {1}" -f $dd.nombre, $_.Exception.Message)
    }
}

# --- 1. El dato manda -------------------------------------------------------
$ult = Get-ChildItem $destino -Recurse -Filter *.parquet -ErrorAction SilentlyContinue |
       Sort-Object LastWriteTime -Descending | Select-Object -First 1
$edad = if ($ult) { ((Get-Date) - $ult.LastWriteTime).TotalMinutes } else { 9999 }

# --- 1.bis. Parte de estado, de un vistazo ---------------------------------
# El fallo de este proyecto nunca ha sido no detectar; ha sido detectar y que
# el aviso se quede en un archivo que nadie abre. Este parte se reescribe en
# CADA pasada, asi que su propia fecha ya dice si el vigilante vive.
try {
    $lin = @("PARTE DEL VIGILANTE -- " + (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), "")
    $lin += "  captura : ultimo volcado hace {0,6:N1} min   {1}" -f $edad,
            $(if ($edad -lt $MIN_SIN_VOLCAR) { "OK" } else { "*** SIN VOLCAR: se repone ***" })
    foreach ($dd in $destinos) {
        $m = Join-Path $dir ("telemetria\" + $dd.marca)
        $txt = "nunca"
        if (Test-Path $m) {
            $tt = [datetime]::MinValue
            if ([datetime]::TryParse((Get-Content $m -Raw).Trim(), [ref]$tt)) {
                $h = ((Get-Date) - $tt).TotalHours
                $txt = "hace {0,5:N1} h" -f $h
            }
        }
        $pres = if (Test-Path $dd.raiz) { "presente" } else { "*** AUSENTE ***" }
        $lin += "  resp. {0,-6}: {1,-14} {2,-9} -> {3}" -f $dd.nombre, $txt, $pres, $dd.destino
    }
    $lin | Out-File -FilePath (Join-Path $dir "telemetria\ESTADO.txt") -Encoding utf8
} catch { }

if ($edad -lt $MIN_SIN_VOLCAR) { exit 0 }   # sano: ni se registra, para no llenar el log

# --- 2. No arrancar dos veces seguidas -------------------------------------
if (Test-Path $marca) {
    $t = Get-Content $marca -Raw
    $tt = [datetime]::MinValue
    if ([datetime]::TryParse($t.Trim(), [ref]$tt)) {
        if (((Get-Date) - $tt).TotalMinutes -lt $MIN_TRAS_ARRANCAR) {
            Registro ("ESPERANDO: arranque hace {0:N1} min, aun sin volcar. No se toca." -f ((Get-Date) - $tt).TotalMinutes)
            exit 0
        }
    }
}

# --- 3. Matar lo que quede colgado, si se puede saber cual es ---------------
# Best effort: si WMI no responde, se sigue igual. Un proceso colgado que no
# vuelca no produce dato, asi que reponerlo es correcto aunque el viejo siga.
try {
    Get-CimInstance Win32_Process -Filter "Name='python.exe'" -ErrorAction Stop |
        Where-Object { $_.CommandLine -like "*captura_estacional*" } |
        ForEach-Object {
            Registro ("COLGADO: se mata PID {0} ({1:N1} min sin volcar)" -f $_.ProcessId, $edad)
            Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
        }
} catch {
    Registro ("AVISO: no se pudo consultar la tabla de procesos ({0}). Se repone igual." -f $_.Exception.GetType().Name)
}

# --- 4. Reponer -------------------------------------------------------------
$log = Join-Path $dir ("telemetria\estacional_{0}.log" -f (Get-Date -Format "yyyyMMdd_HHmmss"))
try {
    $p = Start-Process -FilePath $py `
            -ArgumentList @("captura_estacional.py", "--dias=21", "--sin-guarda-v33") `
            -WorkingDirectory $dir `
            -RedirectStandardOutput $log `
            -RedirectStandardError ($log -replace "\.log$", ".err") `
            -WindowStyle Hidden -PassThru
    (Get-Date).ToString("o") | Out-File -FilePath $marca -Encoding ascii
    Registro ("REPUESTO: PID {0} tras {1:N1} min sin volcar" -f $p.Id, $edad)
} catch {
    Registro ("FALLO AL ARRANCAR: {0}" -f $_.Exception.Message)
    exit 1
}
