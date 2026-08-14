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

# --- 0. Respaldo incremental, MEJOR ESFUERZO --------------------------------
# Va primero y entero dentro de try/catch: si el respaldo falla, el vigilante
# tiene que seguir haciendo su trabajo. Un respaldo caido es una molestia; una
# captura caida son horas de mercado.
#
# `/XO` copia solo lo que no esta ya, y los parquet no se reescriben nunca, asi
# que en regimen son ~2 archivos por pasada. Se estrangula a una vez por hora
# para no releer 1900 nombres de archivo cada 5 min.
#
# [!] EL RESPALDO ESTA EN EL MISMO DISCO FISICO. Protege contra borrado
#     accidental, contra un script mio que la lie y contra una escritura a
#     medias; NO protege contra fallo del SSD. Eso exige un medio aparte y no
#     hay ninguno en esta maquina.
try {
    $marcaResp = Join-Path $dir "telemetria\ultimo_respaldo.txt"
    $toca = $true
    if (Test-Path $marcaResp) {
        $tr = [datetime]::MinValue
        if ([datetime]::TryParse((Get-Content $marcaResp -Raw).Trim(), [ref]$tr)) {
            $toca = ((Get-Date) - $tr).TotalMinutes -ge 60
        }
    }
    if ($toca) {
        $null = robocopy (Join-Path $dir "telemetria") `
                    "C:\Users\Usuario\respaldo_micelio\telemetria" `
                    /E /XO /R:1 /W:2 /NFL /NDL /NP /NJH /NJS /MT:4
        # robocopy devuelve 0-7 en exito (1 = se copio algo). >= 8 es error.
        if ($LASTEXITCODE -ge 8) {
            Registro ("RESPALDO: robocopy devolvio {0}" -f $LASTEXITCODE)
        } else {
            (Get-Date).ToString("o") | Out-File -FilePath $marcaResp -Encoding ascii
        }
    }
} catch {
    Registro ("RESPALDO fallido (no bloquea): {0}" -f $_.Exception.Message)
}

# --- 1. El dato manda -------------------------------------------------------
$ult = Get-ChildItem $destino -Recurse -Filter *.parquet -ErrorAction SilentlyContinue |
       Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($ult) {
    $edad = ((Get-Date) - $ult.LastWriteTime).TotalMinutes
    if ($edad -lt $MIN_SIN_VOLCAR) { exit 0 }   # sano: ni se registra, para no llenar el log
} else {
    $edad = 9999
}

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
