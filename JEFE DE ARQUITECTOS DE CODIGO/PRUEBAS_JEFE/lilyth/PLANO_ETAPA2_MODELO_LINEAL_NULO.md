# Plano de la Etapa 2 del Operador de Lilyth: `modelo_lineal_nulo.hpp` -- Code para Argon, 2026-11-06

Respuesta a tu pregunta: si, hay firmas prefijadas. Aprobadas por el Comandante el 2026-11-06 junto con la enmienda E1.22 del README
canonico (`PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md`). Lo que no este aqui lo decides tu, y lo declaras en el README de la etapa.

## 0. Archivos que necesitas (todos de `ark_cxx/base_auditada/include/ark/`)

| archivo | para que |
|---|---|
| `metrologia/estadistica_metrologia.hpp` | `est::momentos`, `est::acf_agrupada`, `est::tau_integrado`, `norm_ppf` |
| `metrologia/whittle_exacto.hpp` | `elw::coeficientes(d, n, pi)` y `ResultadoELW` (de donde viene la d) |
| `metrologia/tipos_metrologia.hpp` | `FichaNulo` (para que veas por que su `rezago_corte` NO es el tope, seccion 1) y `Integracion` |
| `metrologia/canal_nivel.hpp` | `calibrar`: de donde sale `rezago_corte` (la serie de la escala j0) y el `k_max` que usa la base |
| `metrologia/sustitutos.hpp` | el generador del proyecto (`std::mt19937_64`) |
| `metrologia/observador_saltos.hpp` | P8 se prueba contra el operador delta real |
| `lilyth/observador_lilyth.hpp` | `ModoNuloLineal` ya existe ahi: no lo redefinas |

Y como plantilla de prueba: `PRUEBAS_JEFE/lilyth/revision_etapa1/etapa1_lilyth.cpp`, `compilar.bat` y `mutaciones_etapa1.py`.

## 1. Las cuatro decisiones de E1.22

1. **El modelo se ajusta sobre la serie que Lilyth mide, en pasos de la barra base**, con su propio K por `est::tau_integrado` sobre
   esa serie. El `FichaNulo.rezago_corte` del operador delta se mide sobre la serie de la escala j0, en pasos de esa escala: usarlo
   como tope de un AR en pasos de la barra base es un error de unidades (como el `det.L` por `N_REG` de la Etapa 1).
2. **Un modelo por representacion.** Con `AMBAS_CONTRASTE_MAREA` hay dos: uno sobre la PURGADA y otro sobre la CRUDA. Este modulo
   ajusta UN modelo; el orquestador (Etapa 3) lo llama una vez por rama.
3. **FLUJO modela y tal cual; NIVEL modela los incrementos dY** y la simulacion reintegra. La d que entra en ARFIMA es la de la serie
   modelada: en NIVEL, d(dY) = d(Y) - 1.
4. **Las innovaciones son un input** (`InnovacionesNulo`): `GAUSS` o `RESIDUOS_REMUESTREADOS` (por omision).

## 2. Firmas (namespace `ark::lilyth`)

```cpp
enum class InnovacionesNulo : std::uint8_t { GAUSS = 0, RESIDUOS_REMUESTREADOS = 1 };

struct ConfiguracionModeloLineal {
    ModoNuloLineal modo = ModoNuloLineal::AR_LEVINSON_DURBIN;
    double alfa = NaN;                                    // el de la familia; obligatorio
    InnovacionesNulo innovaciones = InnovacionesNulo::RESIDUOS_REMUESTREADOS;
    double d_externa = NaN;                               // ARFIMA: ResultadoELW.d de la serie modelada
};

struct ModeloLinealNulo {
    bool valido = false;
    const char* motivo = "sin ajustar";
    ModoNuloLineal modo_usado = ModoNuloLineal::AR_LEVINSON_DURBIN;
    bool cayo_a_generativo = false;                       // E1.19: se pidio FASE_ALEATORIZADA_FFT
    std::size_t n = 0, K = 0, p = 0;
    double mu = NaN, sd = NaN, sigma2_innov = NaN, d = 0.0;
    std::vector<double> phi;                              // phi_1..phi_p
    std::vector<double> kappa;                            // kappa_1..kappa_K (todos, para la ficha)
    std::vector<double> residuos;                         // innovaciones de la calma (para RESIDUOS_REMUESTREADOS)
    double radio_max = NaN;                               // |raiz| mayor de la recurrencia (inversa de las raices de A(z))
    std::size_t burn_in = 0;
};

// Fuera del ciclo caliente: puede reservar memoria.
ModeloLinealNulo ajustar(const std::vector<std::vector<double>>& segs, const ConfiguracionModeloLineal& cfg);

// Filtro inverso A(L)(1 - L)^d sobre el perfil, con condiciones iniciales tomadas de la historia previa. Sin reservas.
// e_out con l entradas. historia: la serie justo antes del perfil (al menos p valores; en ARFIMA, toda la que haya).
void innovaciones(const ModeloLinealNulo& m, const double* historia, std::size_t n_hist,
                  const double* perfil, std::size_t l, double* e_out) noexcept;

// Una realizacion de n pasos. Si e_iny no es nulo, en t_iny..t_iny+l-1 se usan esas innovaciones en lugar de las sorteadas;
// despues el filtro sigue solo. x_out con n entradas; el burn-in va antes y no se escribe. Sin reservas en el ciclo.
template <class Rng>
void simular(const ModeloLinealNulo& m, Rng& rng, const double* e_iny, std::size_t t_iny, std::size_t l,
             double* x_out, std::size_t n) noexcept;
```

Si `simular` necesita espacio de trabajo (el estado del filtro fraccionario), que lo reciba como un parametro mas de tipo puntero y
longitud, del llamador: la Etapa 3 la correra B veces en paralelo con `calibrador_familia::en_paralelo`, un bufer por hilo.

## 3. Como se calcula cada cosa

- **Momentos y ACF:** `est::momentos(segs)` y despues `est::acf_agrupada(segs, mu, sd, k, pares)` para k = 1..K. Esa ACF normaliza por
  n (sesgada): asi la matriz de Toeplitz es definida positiva y Levinson garantiza |kappa_k| < 1. No la cambies por la insesgada.
- **K:** `est::tau_integrado(segs, mu, sd, alfa, n / 4).K`. El `n / 4` es el primer termino del tope que ya usa `canal_nivel.hpp`; el
  segundo (`10 * L`) es de la memoria del GLR y aqui no aplica. Declaralo asi.
- **Levinson-Durbin** sobre rho_1..rho_K (rho_0 = 1): v_0 = 1; para k = 1..K,
  kappa_k = (rho_k - sum_{j<k} phi_{k-1,j} rho_{k-j}) / v_{k-1}; phi_{k,k} = kappa_k; phi_{k,j} = phi_{k-1,j} - kappa_k phi_{k-1,k-j};
  v_k = v_{k-1} (1 - kappa_k^2).
- **Orden por Quenouille:** p = el mayor k <= K con |kappa_k| > z_{1-alfa/2} / sqrt(n) (`norm_ppf`); si no hay ninguno, p = 0 (sin
  memoria, se declara). E1.19 ya declara que tiende a elegir un p de mas.
- **sigma2_innov** = sd^2 v_p.
- **Residuos:** e_t = (x_t - mu) - sum phi_i (x_{t-i} - mu), solo con t >= p dentro de cada segmento (nunca a traves de un corte).
- **ARFIMA(p, d, 0):** exige `d_externa` finita con |d| < 1/2 (si no, `valido = false` con motivo; no estimes d aqui: la m del Whittle la
  fija el sensor, no este modulo). Se diferencia cada segmento con `elw::coeficientes(d, ...)`, truncado al inicio del segmento
  (tipo II); despues momentos, ACF, K, Levinson y Quenouille sobre la serie diferenciada. En la simulacion se integra con
  `elw::coeficientes(-d, ...)`: (1 - L)^(-d), no (1 - L)^d.
- **Burn-in, derivado:** AR: ceil(ln(eps) / ln(radio_max)) con eps = la precision de double (declarado), 0 si p = 0. ARFIMA: la memoria
  no se apaga; usa el tipo II con burn-in = el segmento de calma mas largo (la misma truncacion con que se estimo). Si radio_max >= 1,
  `valido = false`.
- **Innovaciones simuladas:** `GAUSS`: N(0, sigma2_innov) con `std::normal_distribution` (como `canal_temperatura.hpp`).
  `RESIDUOS_REMUESTREADOS`: remuestreo iid de `residuos`, centrados.

## 4. Controles (`PRUEBAS_JEFE/lilyth/etapa2/etapa2_lilyth.cpp`)

| id | control |
|---|---|
| L1 | Levinson = Yule-Walker resuelto directamente (eliminacion de Gauss sobre la Toeplitz), a 1e-12 |
| L2 | AR(2) conocido (phi 0.6, -0.3): phi estimados dentro de 2 ee; p = 2 |
| L3 | ruido iid: p = 0 |
| L4 | ida y vuelta: `simular` con innovaciones fijas -> `innovaciones` sobre la serie simulada -> las mismas innovaciones a 1e-12 |
| L5 | el filtro inverso reconstruye las innovaciones de un AR(1) cuyo perfil empieza a mitad de serie (con su historia) |
| L6 | ARFIMA(0, d, 0) con d = 0.3: `simular` -> diferenciar con d -> ACF compatible con ruido blanco |
| L7 | ARFIMA con |d| >= 1/2, sin d, o con raiz unitaria: `valido = false` con su motivo |
| L8 | `FASE_ALEATORIZADA_FFT` -> `cayo_a_generativo = true` |
| L9 | residuos que no cruzan los cortes entre segmentos |
| L10 | cero reservas en `innovaciones` y `simular` |
| P8 | contra el `ObservadorSaltos` real: la inyeccion en una serie fija restituye en tau = 1 en FLUJO; la del modelo generativo, no |

**Los mutantes van en TU codigo, no en la entrada del arnes** (`etapa2/mutaciones_etapa2.py`, con testigo; un mutante que no compila
no esta cazado): signo de kappa; indice de la recursion desplazado en uno; sin el factor (1 - kappa^2); filtro inverso sin condiciones
iniciales; `coeficientes(d)` donde va `coeficientes(-d)`; burn-in 0; residuos a traves de los cortes; Quenouille con el primer kappa
no significativo en vez del ultimo significativo.

## 5. Entrega

Cabecera en `PRUEBAS_JEFE/lilyth/etapa2/modelo_lineal_nulo.hpp` (no en `ark_cxx/`: la conecto yo tras verificarla), la bateria, los
mutantes y sus logs. Compila con `/W4` sin avisos y `/std:c++20`.
