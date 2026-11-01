"""
enumeration.py
Enumeracion exacta para SVP y CVP (Fincke-Pohst / Kannan), aritmetica racional.

Fuentes:
  [Pohst81]   M. Pohst, SIGSAM Bulletin 15 (1981), 37-44.
  [FP85]      U. Fincke, M. Pohst, Math. Comp. 44(170) (1985), 463-471. Idea: recorrer los
              enteros x dentro de un ELIPSOIDE (no de una caja).
  [Kannan87]  R. Kannan, "Minkowski's convex body theorem and integer programming",
              Math. Oper. Res. 12(3) (1987).
  [Galbraith] cap. 18.4 y 18.4.1. Usada para la implementación:
      Con v = sum x_i b_i, B_i = ||b*_i||^2, mu_{j,i} = <b_j,b*_i>/<b*_i,b*_i> y
      z_i = x_i + sum_{j>i} mu_{j,i} x_j, el vector cumple ||v||^2 <= A si y solo si
          sum_i z_i^2 B_i <= A,
      y las sumas parciales dan las cotas de cada coordenada (Teorema 18.4.3):
          z_i^2 B_i <= A - sum_{j>i} z_j^2 B_j.
      Variante CVP (Lema 18.4.8): se sustituye z_i por z_i - y_i, con t = sum y_i b*_i.
  Orden de candidatos: por |x_i - c_i| creciente (zig-zag), idea de Schnorr-Euchner que
  [Galbraith 18.4] menciona;

Que garantiza:  la busqueda es completa dentro del radio, asi que si termina sin agotar el
presupuesto el resultado es el optimo (certificado). Coste exponencial en la dimension.
Por eso se parte siempre de una base reducida (LLL): cotas mas estrechas, menos nodos.
"""
import itertools
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import List, Optional, Sequence

from babai import babai_nearest_plane, dist2, round_frac
from lll import Matrix, dot, gso_full, norm2


class BudgetExceeded(Exception):
    pass


@dataclass
class EnumResult:
    vector: Optional[List[int]]   # mejor vector con coste estrictamente menor que el radio (o none)
    coeffs: Optional[List[int]]   # sus coordenadas respecto a la base recibida
    dist2: Optional[int]
    certified: bool               # True <=> la busqueda dentro del radio se completo
    nodes: int
    seconds: float


def _zigzag(x0: int, s: int):
    """x0, x0+s, x0-s, x0+2s, x0-2s, ...
    Propiedad: si x0 = round(c) y d = |c - x0| <= 1/2,
    s = sign(c - x0), las distancias a c son d, 1-d, 1+d, 2-d, 2+d, ... y como d <= 1/2
    cumplen d <= 1-d <= 1+d <= 2-d <= ..."""
    yield x0
    m = 1
    while True:
        yield x0 + m * s
        yield x0 - m * s
        m += 1


def enumerate_vectors(
    B: Matrix,
    target: Optional[Sequence[int]] = None,
    radius2=None,
    max_nodes: Optional[int] = None,
    max_seconds: Optional[float] = None,
) -> EnumResult:
    """
    target=None  -> SVP: vector no nulo de L(B) mas corto.
    target=t     -> CVP: punto de L(B) mas cercano a t.
    radius2      -> solo se buscan vectores con distancia^2 menor.
                    Por defecto: SVP = ||fila mas corta||^2 ; CVP = distancia^2 de Babai.
    Si no existe ninguno mas corto que el radio, devuelve vector=None (el candidato inicial
    ya era optimo) y certified=True si la busqueda fue completa.
    """
    t0 = time.perf_counter()
    n = len(B)
    svp = target is None
    bstar, mu, Bn = gso_full(B)

    if radius2 is None:
        if svp:
            radius2 = min(norm2(r) for r in B)
        else:
            radius2 = dist2(babai_nearest_plane(B, target)[0], target)
    bound2 = Fraction(radius2)

    tau = [Fraction(0)] * n if svp else [dot(target, bstar[i]) / Bn[i] for i in range(n)]
    x = [0] * n
    best_x: Optional[List[int]] = None
    nodes = 0

    def rec(i: int, rho: Fraction, zero_above: bool):
        nonlocal nodes, bound2, best_x
        c = tau[i]
        for j in range(i + 1, n):
            if x[j]:
                c -= mu[j][i] * x[j]
        Bi = Bn[i]
        if svp and zero_above:
            # Simetria v <-> -v. invariante: mirando desde el indice mas alto, el primer
            # coeficiente no nulo es positivo (los de arriba valen 0 => c = 0 => 0,1,2,...)
            cands = itertools.count(0)
        else:
            x0 = round_frac(c)
            cands = _zigzag(x0, 1 if c - x0 >= 0 else -1)
        for xi in cands:
            diff = xi - c
            cost = rho + Bi * diff * diff
            if cost >= bound2:
                break                              # los siguientes candidatos son peores
            nodes += 1
            if max_nodes is not None and nodes > max_nodes:
                raise BudgetExceeded
            if max_seconds is not None and nodes % 256 == 0 and time.perf_counter() - t0 > max_seconds:
                raise BudgetExceeded
            x[i] = xi
            if i == 0:
                if svp and zero_above and xi == 0:
                    continue                       # excluye el vector nulo
                bound2 = cost                      
                best_x = x[:]
            else:
                rec(i - 1, cost, zero_above and xi == 0)
        x[i] = 0

    certified = True
    try:
        rec(n - 1, Fraction(0), True)
    except BudgetExceeded:
        certified = False

    vec = d2 = None
    if best_x is not None:
        vec = [sum(best_x[i] * B[i][j] for i in range(n)) for j in range(len(B[0]))]
        d2 = norm2(vec) if svp else dist2(vec, target)
        assert d2 == bound2, "incoherencia interna: coste acumulado != distancia recalculada"
    return EnumResult(vec, best_x, d2, certified, nodes, time.perf_counter() - t0)


@dataclass
class SolveResult:
    vector: List[int]
    dist2: int
    certified: bool
    nodes: int
    seconds: float
    source: str          # de donde sale la respuesta: 'enumeracion' | 'base reducida' | 'Babai'
    initial_dist2: int   # coste del candidato de partida (fila mas corta o punto de Babai)


def solve_svp(B: Matrix, max_nodes=None, max_seconds=None) -> SolveResult:
    """SVP sobre una base (LLL-reducida): fila mas corta + enumeracion."""
    cand = min(B, key=norm2)
    r = enumerate_vectors(B, None, norm2(cand), max_nodes, max_seconds)
    if r.vector is not None:
        return SolveResult(r.vector, r.dist2, r.certified, r.nodes, r.seconds, "enumeracion", norm2(cand))
    return SolveResult(list(cand), norm2(cand), r.certified, r.nodes, r.seconds, "base reducida", norm2(cand))


def solve_cvp(B: Matrix, t: Sequence[int], max_nodes=None, max_seconds=None) -> SolveResult:
    """CVP: Babai (plano mas cercano) da el radio inicial; la enumeracion lo mejora o lo certifica."""
    p0, _ = babai_nearest_plane(B, t)
    d0 = dist2(p0, t)
    r = enumerate_vectors(B, t, d0, max_nodes, max_seconds)
    if r.vector is not None:
        return SolveResult(r.vector, r.dist2, r.certified, r.nodes, r.seconds, "enumeracion", d0)
    return SolveResult(p0, d0, r.certified, r.nodes, r.seconds, "Babai", d0)
