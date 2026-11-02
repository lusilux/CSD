# Solucionador de celosías, LWE y Ring-LWE

Implementación propia (en Python, con fracciones exactas) para resolver problemas de celosías (SVP, CVP), LWE y Ring-LWE mediante el pipeline **LLL → Babai → enumeración**. Cada resultado se verifica de forma independiente.

Trabajo de Criptología Post Cuántica, Máster Universitario en Ciberseguridad y Ciberinteligencia (UPV).

## Punto de entrada

`solve.py` es el programa principal. Identifica el tipo de instancia, selecciona el solucionador, mide tiempos y exporta los resultados a CSV.

```
solve.py
 ├─ SVP/CVP → instance_io.py
 │            → lll.py
 │            → babai.py (solo CVP)
 │            → enumeration.py
 │            → verify.py
 │
 ├─ LWE     → lwe_io.py
 │            → lwe_solve.py
 │            → lwe_attack.py
 │            → lll.py + enumeration.py
 │
 └─ RLWE    → rlwe.py
              → lwe_attack.py
              → enumeration.py
```

## Módulos de algoritmos

| Script | Función |
|---|---|
| `lll.py` | Reducción LLL con fracciones exactas (δ = 0.99), determinantes y operaciones de celosías. Comprueba que la base resultante es LLL-reducida y genera la misma celosía. |
| `babai.py` | Algoritmo del plano más cercano de Babai: solución aproximada de CVP, usada para fijar la cota inicial de la enumeración. |
| `enumeration.py` | Enumeración exacta Fincke–Pohst/Kannan con orden de visita de Schnorr–Euchner, para SVP y CVP. Si termina dentro del presupuesto, el resultado está certificado como óptimo. |

## Lectura de instancias

| Script | Función |
|---|---|
| `instance_io.py` | Lectura de instancias SVP/CVP (formato PoliformaT). |
| `lwe_io.py` | Lectura de instancias LWE (formato PoliformaT). |

## Ataques LWE y Ring-LWE

| Script | Función |
|---|---|
| `lwe_attack.py` | Algoritmo matemático del ataque LWE: construye la celosía q-aria, aplica la incrustación de Kannan, recupera el secreto `s` y lo contrasta con una resolución directa de CVP. |
| `lwe_solve.py` | Orquesta el ataque LWE y da formato a la salida. |
| `rlwe.py` | Ataque a Ring-LWE: matriz negacíclica, celosía de dimensión 2n+1, recuperación de `(s, e)` y verificación por multiplicación polinómica independiente (y descifrado de bits si hay criptograma). |

## Verificación

`verify.py` se usa en la ejecución principal para **SVP/CVP**. Comprueba que:

- el vector pertenece a la celosía;
- el vector SVP no es nulo;
- se cumple la desigualdad de Hermite para SVP;
- las coordenadas respecto a la base original son correctas.

En **LWE y Ring-LWE** la verificación está implementada dentro de `lwe_attack.py` y `rlwe.py` (el secreto recuperado debe cumplir las ecuaciones originales y descifrar los bits del criptograma), no mediante `verify.py`.

## Pruebas, tablas y baterías de benchmark

| Script | Función |
|---|---|
| `test_lwe.py` | Pruebas independientes del ataque LWE. |
| `make_tables.py` | Genera las tablas de resultados para la memoria. |
| `bench_run.py` | Batería de problemas aleatorios SVP/CVP (familias q-aria y uniforme) de dimensión creciente. |
| `lwe_bench.py` | Batería de problemas LWE aleatorios de dimensión creciente. |
| `plot_bench.py` | Genera las figuras a partir de los resultados de las baterías. |

## Uso

<!-- Completar con los argumentos reales de solve.py, por ejemplo: -->

```bash
python solve.py <ficheros_de_instancia> [opciones]
```

Parámetros de entrada de `solve.py`: ficheros de instancia, δ, backend y presupuesto (nodos, segundos). La salida incluye respuesta, estado (`OK_CERTIFICADO`, `SIN_CERTIFICAR`, `FALLO_VERIFICACION`), métricas por instancia y un fichero CSV.
