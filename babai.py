"""
babai.py -- Algoritmos de Babai para CVP aproximado, con aritmetica exacta.

Fuentes:
  [Babai86]   L. Babai, "On Lovasz' lattice reduction and the nearest lattice point
              problem", Combinatorica 6 (1986), 1-13.
  [Galbraith] S. Galbraith, "Mathematics of Public Key Cryptography", cap. 18:
              18.1 (Algoritmo 26, plano mas cercano) y 18.2 (redondeo).
  [PDF]       Enunciado de la asignatura, seccion GGH (descifrado): y = Round(c B^-1).

Convenio: bases por FILAS; el resultado es un punto de la celosia (lista de enteros) y sus
coordenadas enteras respecto a la base recibida.

Garantias (Galbraith, con base LLL-reducida y delta=3/4):
  plano mas cercano:  ||v - t|| < 2^(n/2) * dist(t, L)      (Teorema 18.1.6)
  redondeo:           dependiente de la base, peor cota (Teorema 18.2.1)
  Ninguno garantiza el punto MAS cercano; para eso esta la enumeracion.
"""
import math
from fractions import Fraction
from typing import List, Sequence, Tuple

from lll import Matrix, dot, gso_full, norm2, solve_row


def round_frac(x: Fraction) -> int:
    """Entero mas cercano (empates hacia arriba)."""
    return math.floor(x + Fraction(1, 2))


def _combine(B: Matrix, x: Sequence[int]) -> List[int]:
    return [sum(x[i] * B[i][j] for i in range(len(B))) for j in range(len(B[0]))]


def dist2(u: Sequence, v: Sequence) -> int:
    return norm2([a - b for a, b in zip(u, v)])


def babai_round(B: Matrix, t: Sequence[int]) -> Tuple[List[int], List[int]]:
    """
    Redondeo de Babai [Galbraith 18.2]: t = sum l_i b_i  (l = t B^-1, exacto),
    salida v = sum round(l_i) b_i.  Es el 'Babai round' del descifrado GGH del PDF.
    Devuelve (v, x) con v = x B.
    """
    l = solve_row(B, t)
    x = [round_frac(c) for c in l]
    return _combine(B, x), x


def babai_nearest_plane(B: Matrix, t: Sequence[int]) -> Tuple[List[int], List[int]]:
    """
    Plano mas cercano de Babai [Babai86; Galbraith 18.1, Alg. 26].
    Para i = n..1:  l_i = <w, b*_i> / <b*_i, b*_i>,  x_i = round(l_i),  w <- w - x_i b_i.
    (Galbraith resta ademas (l_i - x_i) b*_i; esa componente es ortogonal a todos los
    b*_j con j < i, asi que no cambia los l_j posteriores y es equivalente.)
    Devuelve (v, x) con v = x B.
    """
    bstar, _, Bn = gso_full(B)
    n = len(B)
    w = [Fraction(a) for a in t]
    x = [0] * n
    for i in range(n - 1, -1, -1):
        xi = round_frac(dot(w, bstar[i]) / Bn[i])
        x[i] = xi
        w = [wj - xi * B[i][j] for j, wj in enumerate(w)]
    return _combine(B, x), x
