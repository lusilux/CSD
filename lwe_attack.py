"""
lwe_attack.py -- Ataque a LWE (Regev) por reduccion a celosias: celosia q-aria + incrustacion
de Kannan, resuelto con los MISMOS bloques que SVP/CVP (LLL exacto + enumeracion).

Modelo (PDF, seccion LWE): b_i = <a_i, s> + e_i (mod q),  i = 1..m,  con errores e_i pequenos.
   cifrado de un bit con la clave publica: (u, v) con u = suma de ecuaciones elegidas,
   v = suma de sus b_i + bit*floor(q/2)   (mod q).
   descifrado: d = v - <u, s> (mod q), centrado:  |d| pequeno -> bit 0 ; |d| ~ q/2 -> bit 1.

Algoritmo (cada paso con su fuente; las referencias las doy de memoria, VERIFICAR antes de citarlas):
  1. Celosia q-aria  L = { y in Z^m : y = A s (mod q) para algun s }  (dimension m, det = q^(m-k)).
       Si las k primeras filas A1 son invertibles mod q:  y2 = (A2 A1^-1) y1 (mod q), con base
           [ I_k | (A2 A1^-1)^T ]    y    [ 0 | q I_(m-k) ]
       Es la construccion estandar de celosia q-aria [Micciancio-Regev 2009, "Lattice-based
       cryptography"; Regev 2005/2009].
  2. b esta a distancia ||e|| de L: resolver LWE == resolver CVP con error pequeno (BDD) en L.
  3. Incrustacion de Kannan [Kannan 1987] (CVP -> SVP): a la base de L (con una columna extra de
       ceros) se le anade la fila (b, M).  El vector (b - y, M) = (e, M) esta en la celosia
       incrustada y es MUY corto -> SVP unico.  Con M = 1 es la eleccion habitual en el ataque
       'primal' [Albrecht-Player-Scott 2015, "On the concrete hardness of LWE"].
  4. LLL + enumeracion (modulos lll.py y enumeration.py, ya existentes) -> vector mas corto.
  5. De e se recupera s resolviendo A1 s = b1 - e1 (mod q); se comprueban TODAS las ecuaciones.
  6. Verificacion cruzada: CVP directo (Babai + enumeracion) en L sin incrustar; debe dar el mismo s.

Que se certifica: la enumeracion completa demuestra que (e, M) es el vector no nulo mas corto de la
celosia incrustada; entonces ningun otro punto de L esta mas cerca de b que y (cualquier y' mas
cercano daria un vector (b-y', M) mas corto), es decir, e es el error de norma minima.
NO se demuestra que coincida con el error 'real' del cifrador (eso solo es cierto con alta
probabilidad cuando el error es pequeno frente a lambda_1(L)); se mide con el 'hueco' gauss/||(e,M)||.
"""
import math
import time
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

from enumeration import solve_cvp, solve_svp
from lll import Matrix, Trace, det_int, is_lll_reduced, lll_reduce, norm2, same_lattice, same_lattice_det

EXACT_MAX_DIM_LWE = 40   # por encima se usa fpylll si esta instalado (umbral provisional)


# ------------------------------------------------------------------ algebra modular
def centered(a: int, q: int) -> int:
    a %= q
    return a - q if a > q // 2 else a


def _unit(a: int, q: int) -> bool:
    return math.gcd(a % q, q) == 1


def pick_independent_rows(A: Matrix, q: int) -> List[int]:
    """Indices de k filas de A cuya submatriz k x k es invertible mod q (eliminacion modular)."""
    k = len(A[0])
    basis: List[Tuple[int, List[int]]] = []          # (columna pivote, fila reducida con pivote 1)
    chosen: List[int] = []
    for i, row in enumerate(A):
        r = [x % q for x in row]
        for c, pr in basis:
            if r[c]:
                f = r[c]
                r = [(x - f * y) % q for x, y in zip(r, pr)]
        c = next((j for j in range(k) if _unit(r[j], q)), None)
        if c is None:
            continue
        inv = pow(r[c], -1, q)
        r = [(x * inv) % q for x in r]
        basis.append((c, r))
        chosen.append(i)
        if len(chosen) == k:
            return chosen
    raise ValueError("las ecuaciones no tienen rango k modulo q (o q no es primo y falta una unidad)")


def inv_mod_matrix(M: Matrix, q: int) -> Matrix:
    n = len(M)
    A = [[x % q for x in row] + [int(i == j) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next((i for i in range(c, n) if _unit(A[i][c], q)), None)
        if p is None:
            raise ValueError("matriz no invertible mod q")
        A[c], A[p] = A[p], A[c]
        inv = pow(A[c][c], -1, q)
        A[c] = [(x * inv) % q for x in A[c]]
        for i in range(n):
            if i != c and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % q for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


def mat_vec_mod(M: Matrix, v: Sequence[int], q: int) -> List[int]:
    return [sum(a * b for a, b in zip(row, v)) % q for row in M]


def mat_mul_mod(X: Matrix, Y: Matrix, q: int) -> Matrix:
    return [[sum(X[i][t] * Y[t][j] for t in range(len(Y))) % q for j in range(len(Y[0]))]
            for i in range(len(X))]


# ------------------------------------------------------------------ celosia q-aria
def qary_lattice_basis(A: Matrix, q: int):
    """
    Base (por filas, m x m) de L = {A s mod q} + q Z^m, con las filas de A reordenadas para que
    las k primeras sean independientes. Devuelve (base, perm) con perm[j] = indice original de la
    fila que ocupa la posicion j.
    """
    m, k = len(A), len(A[0])
    P = pick_independent_rows(A, q)
    perm = P + [i for i in range(m) if i not in set(P)]
    Ap = [A[i] for i in perm]
    A1inv = inv_mod_matrix(Ap[:k], q)
    Mmat = mat_mul_mod(Ap[k:], A1inv, q) if m > k else []          # (m-k) x k
    basis: Matrix = []
    for i in range(k):                                              # (e_i | M[:, i])
        basis.append([int(i == j) for j in range(k)] + [Mmat[r][i] for r in range(m - k)])
    for j in range(m - k):                                          # (0 | q e_j)
        basis.append([0] * k + [q if r == j else 0 for r in range(m - k)])
    assert abs(det_int(basis)) == q ** (m - k)
    return basis, perm


def _reduce(B: Matrix, delta: Fraction, backend: str, stats: Optional[dict] = None):
    """LLL + comprobaciones. Si se pasa `stats` (dict) se rellena con t_lll, t_check y swaps
    (swaps solo con el backend exacto; fpylll no los expone -> None)."""
    n = len(B)
    if backend == "auto":
        backend = "exact" if n <= EXACT_MAX_DIM_LWE else "fpylll"
    if backend == "exact":
        tr = Trace(False)                                   # solo cuenta intercambios
        t = time.perf_counter()
        Bred, U = lll_reduce(B, delta, tr)
        t_lll = time.perf_counter() - t
        t = time.perf_counter()
        assert is_lll_reduced(Bred, delta) and same_lattice(B, U, Bred)
        t_check = time.perf_counter() - t
        if stats is not None:
            stats.update(t_lll=t_lll, t_check=t_check, swaps=tr.swaps)
        return Bred, "exact"
    from fpylll import LLL, IntegerMatrix            # importacion perezosa
    t = time.perf_counter()
    A = IntegerMatrix.from_matrix(B)
    LLL.reduction(A, delta=float(delta), eta=0.501)
    Bred = [[A[i, j] for j in range(A.ncols)] for i in range(A.nrows)]
    t_lll = time.perf_counter() - t
    t = time.perf_counter()
    assert same_lattice_det(B, Bred)
    if stats is not None:
        stats.update(t_lll=t_lll, t_check=time.perf_counter() - t, swaps=None)
    return Bred, "fpylll"


# ------------------------------------------------------------------ ataque
def solve_secret(A: Matrix, b: Sequence[int], q: int, M: int = 1, delta: Fraction = Fraction(99, 100),
                 backend: str = "auto", samples: Optional[int] = None, cross_check: bool = True,
                 max_nodes=None, max_seconds=None) -> Dict:
    """Recupera (s, e) con b = A s + e (mod q), e de norma minima. Ver docstring del modulo."""
    m_all, k = len(A), len(A[0])
    m = samples or m_all
    if not (k <= m <= m_all):
        raise ValueError(f"samples debe estar entre k={k} y m={m_all}")
    A, b = A[:m], list(b[:m])
    out = dict(ok=False, note="", m_used=m, k=k, q=q, dim=m + 1)

    t0 = time.perf_counter()
    L, perm = qary_lattice_basis(A, q)
    bp = [b[i] % q for i in perm]
    E = [row + [0] for row in L] + [bp + [M]]                    # incrustacion de Kannan
    st: dict = {}
    Ered, used = _reduce(E, delta, backend, st)
    out["t_reduce"], out["backend"] = time.perf_counter() - t0, used
    out.update(t_lll=st["t_lll"], t_check=st["t_check"], swaps=st["swaps"])

    def direct_cvp():                                            # CVP exacto en L (sin incrustar)
        Lred, _ = _reduce(L, delta, backend)
        c = solve_cvp(Lred, bp, max_nodes, max_seconds)
        return c, [bp[j] - c.vector[j] for j in range(m)]

    r = solve_svp(Ered, max_nodes, max_seconds)
    w = list(r.vector)
    if w[-1] < 0:
        w = [-x for x in w]
    out.update(t_solve=r.seconds, nodes=r.nodes, certified=r.certified, source=r.source,
               value=norm2(w), method="incrustacion de Kannan + SVP")
    cvp_cache = None
    if w[-1] != M:
        # L tiene un vector mas corto que (e, M): la incrustacion no sirve, pero el CVP directo
        # (Babai + enumeracion) sigue siendo exacto. Se usa como respaldo.
        t1 = time.perf_counter()
        c, e_perm = direct_cvp()
        out.update(method="CVP directo (respaldo: la incrustacion no aisla el error)", certified=c.certified,
                   t_solve=out["t_solve"] + time.perf_counter() - t1, nodes=out["nodes"] + c.nodes,
                   source=c.source, value=norm2(e_perm) + M * M)
        cvp_cache = (c, e_perm)
    else:
        e_perm = w[:m]
    e = [0] * m
    for j, i in enumerate(perm):
        e[i] = e_perm[j]
    P = perm[:k]
    A1inv = inv_mod_matrix([A[i] for i in P], q)
    s = mat_vec_mod(A1inv, [(b[i] - e[i]) % q for i in P], q)
    res = [centered(b[i] - sum(a * x for a, x in zip(A[i], s)), q) for i in range(m)]
    consistent = all((res[i] - e[i]) % q == 0 for i in range(m))
    out.update(secret=s, errors=e, errors_centered=[centered(x, q) for x in e], consistent=consistent,
               err_norm2=norm2(e), err_max=max(abs(x) for x in e))
    # hueco: lambda_1(L) estimado (heuristica gaussiana) frente al error
    gh = math.sqrt(m / (2 * math.pi * math.e)) * float(q) ** ((m - k) / m)
    out["gh_lambda1"] = gh
    out["gap"] = gh / math.sqrt(norm2(e_perm)) if norm2(e_perm) else float("inf")
    # soluciones espurias esperadas: s' en Z_q^k con TODOS los errores en [-B, B]  (A aleatoria):
    #   q^k * ((2B+1)/q)^m ; si es << 1 la solucion es unica con probabilidad ~ 1
    B_ = out["err_max"]
    out["spurious_log10"] = k * math.log10(q) + m * math.log10(min(1.0, (2 * B_ + 1) / q))

    if cvp_cache is not None:
        out["cvp_agrees"] = True                                  # el resultado YA es el del CVP directo
        out["cvp_err_norm2"] = norm2(e_perm)
    elif cross_check:                                             # verificacion cruzada con CVP directo
        c, e2 = direct_cvp()
        out["cvp_agrees"] = (c.certified and norm2(e2) == norm2(e_perm) and
                             [centered(x, q) for x in e2] == [centered(x, q) for x in e_perm])
        out["cvp_err_norm2"] = norm2(e2)
    out["ok"] = consistent and out["certified"]
    return out


def decrypt_bit(u: Sequence[int], v: int, s: Sequence[int], q: int) -> Dict:
    R = sum(a * x for a, x in zip(u, s)) % q
    d = centered(v - R, q)                      # ~ error acumulado (bit 0) o ~ q/2 + error (bit 1)
    half = q // 2
    d1 = centered(v - R - half, q)              # distancia al polo del bit 1
    return dict(R=R, d=d, bit=0 if abs(d) < abs(d1) else 1, dist0=abs(d), dist1=abs(d1), half=half)
