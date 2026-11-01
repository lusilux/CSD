"""
verify.py -- Verificacion INDEPENDIENTE de las respuestas (no reutiliza la enumeracion ni Babai).

  - pertenencia a la celosia (coordenadas enteras respecto a la base ORIGINAL),
  - desigualdad de Hermite [Nguyen, Teorema 8]: lambda1^2 <= (4/3)^((n-1)/2) det^(2/n),
    condicion NECESARIA para que un vector sea el mas corto (sirve para cazar errores gordos),
  - busqueda por fuerza bruta en una caja de coeficientes (solo dimension pequena; NO es una
    prueba, es un contraste con un algoritmo completamente distinto).
"""
import itertools
import math
from fractions import Fraction
from typing import List, Optional, Sequence

from babai import babai_nearest_plane, dist2
from lll import Matrix, det_int, in_lattice, norm2, solve_row


def coefficients(B: Matrix, v: Sequence) -> Optional[List[int]]:
    """Coordenadas enteras de v respecto a B (o None si v no esta en la celosia)."""
    x = solve_row(B, v)
    if x is None or any(c.denominator != 1 for c in x):
        return None
    return [int(c) for c in x]


def hermite_ok(B: Matrix, v: Sequence) -> bool:
    """||v||^(2n) <= (4/3)^(n(n-1)/2) * det^2   (exacto, en enteros/fracciones)."""
    n = len(B)
    d = det_int(B)
    return Fraction(norm2(v)) ** n <= Fraction(4, 3) ** (n * (n - 1) // 2) * d * d


def gaussian_heuristic2(B: Matrix) -> float:
    """Estimacion (heuristica, NO una cota) de lambda1^2: (n/(2 pi e))^1 * det^(2/n)."""
    n = len(B)
    return n / (2 * math.pi * math.e) * abs(det_int(B)) ** (2.0 / n)


def verify_svp(B: Matrix, v: Sequence) -> dict:
    ok_lat = coefficients(B, v) is not None
    return dict(in_lattice=ok_lat, nonzero=any(v), hermite=hermite_ok(B, v),
                ok=ok_lat and any(v) and hermite_ok(B, v))


def verify_cvp(B: Matrix, t: Sequence, p: Sequence) -> dict:
    ok_lat = coefficients(B, p) is not None
    return dict(in_lattice=ok_lat, ok=ok_lat)


def brute_force_svp(B: Matrix, K: int):
    """Mejor vector no nulo con coeficientes en [-K, K] respecto a B (B conviene reducida)."""
    n, m = len(B), len(B[0])
    best = None
    for x in itertools.product(range(-K, K + 1), repeat=n):
        if not any(x):
            continue
        v = [sum(x[i] * B[i][j] for i in range(n)) for j in range(m)]
        d = norm2(v)
        if best is None or d < best[0]:
            best = (d, v)
    return best


def brute_force_cvp(B: Matrix, t: Sequence, K: int):
    """Mejor punto con coeficientes en [x0-K, x0+K], x0 = coordenadas del punto de Babai."""
    n, m = len(B), len(B[0])
    _, x0 = babai_nearest_plane(B, t)
    best = None
    for dx in itertools.product(range(-K, K + 1), repeat=n):
        x = [a + b for a, b in zip(x0, dx)]
        v = [sum(x[i] * B[i][j] for i in range(n)) for j in range(m)]
        d = dist2(v, t)
        if best is None or d < best[0]:
            best = (d, v)
    return best
