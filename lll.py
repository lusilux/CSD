"""
lll.py
LLL con aritmetica exacta (con fracciones), con traza y verificadores.

Fuentes usadas (cada bloque indica la suya):
  [LLL82]  A.K. Lenstra, H.W. Lenstra Jr., L. Lovasz, "Factoring polynomials with
           rational coefficients", Math. Ann. 261 (1982), 515-534.
           Algoritmo de reduccion de bases: pseudocodigo reproducido en la fig. 1.13
           de [Smeets]  (version iterativa con actualizacion incremental de la GSO).
  [Nguyen] P.Q. Nguyen, "Hermite's constant and lattice algorithms", cap. 2 de
           Nguyen & Vallee (eds.), The LLL Algorithm, Springer 2010.
           Definiciones 12 (GSO), 15 (size-reduced), LLL-reducida (factor delta),
           Teorema 9 (garantias), Algoritmo 3 (size-reduction), Algoritmo 6 (LLL basico).
  [Smeets] I. Smeets et al., "The history of the LLL-algorithm", cap. 1 del mismo libro.
  [Cohen]  H. Cohen, "A course in computational algebraic number theory", Alg. 2.6.3
           (formulas de actualizacion de mu tras un intercambio; ver nota en swap).
"""
from fractions import Fraction
from typing import List, Optional, Sequence, Tuple

Matrix = List[List[int]]

# utilidades
def dot(u: Sequence, v: Sequence):
    return sum(a * b for a, b in zip(u, v))


def norm2(v: Sequence):
    return dot(v, v)


def det_int(M: Matrix) -> int:
    """Determinante exacto. Método Bareiss"""
    n = len(M)
    A = [row[:] for row in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def gso_full(B: Matrix) -> Tuple[List[List[Fraction]], List[List[Fraction]], List[Fraction]]:
    """
    Ortogonalizacion de Gram-Schmidt exacta  [Nguyen, Def. 12, ec. (2.9)]:
        b*_i = b_i - sum_{j<i} mu_ij b*_j ,   mu_ij = <b_i, b*_j> / ||b*_j||^2
    Devuelve (bstar, mu, Bn): los vectores b*_i, mu[i][j] (j<i) y Bn[i] = ||b*_i||^2.
    """
    n = len(B)
    bstar: List[List[Fraction]] = []
    mu = [[Fraction(0)] * n for _ in range(n)]
    Bn: List[Fraction] = []
    for i in range(n):
        v = [Fraction(x) for x in B[i]]
        for j in range(i):
            mu[i][j] = dot(B[i], bstar[j]) / Bn[j]
            v = [vi - mu[i][j] * bj for vi, bj in zip(v, bstar[j])]
        bstar.append(v)
        Bn.append(norm2(v))
        if Bn[i] == 0:
            raise ValueError("los vectores de entrada no son linealmente independientes")
    return bstar, mu, Bn


def gso(B: Matrix) -> Tuple[List[List[Fraction]], List[Fraction]]:
    """(mu, Bn) de gso_full; se mantiene por compatibilidad."""
    _, mu, Bn = gso_full(B)
    return mu, Bn


# verificadores independientes del algoritmo (recalculan la gso desde cero)
def is_size_reduced(B: Matrix) -> bool:
    """[Nguyen, Def. 15, ec. (2.19)]: |mu_ij| <= 1/2 para todo j < i."""
    mu, _ = gso(B)
    n = len(B)
    return all(abs(mu[i][j]) <= Fraction(1, 2) for i in range(n) for j in range(i))


def is_lll_reduced(B: Matrix, delta: Fraction = Fraction(3, 4)) -> bool:
    """
    [Nguyen, seccion 'The LLL algorithm']: size-reduced y condicion de Lovasz
        ||b*_{i+1} + mu_{i+1,i} b*_i||^2 >= delta ||b*_i||^2
    que equivale a  ||b*_{i+1}||^2 >= (delta - mu_{i+1,i}^2) ||b*_i||^2.
    """
    mu, Bn = gso(B)
    n = len(B)
    if not all(abs(mu[i][j]) <= Fraction(1, 2) for i in range(n) for j in range(i)):
        return False
    return all(Bn[i + 1] >= (delta - mu[i + 1][i] ** 2) * Bn[i] for i in range(n - 1))


def same_lattice(B: Matrix, U: Matrix, Bred: Matrix) -> bool:
    """
    Comprueba que Bred = U * B con U entera y det(U) = +-1 (cambio de base
    unimodular, [Nguyen, Teorema 3]). Bred genera exactamente la misma celosia.
    """
    n, m = len(B), len(B[0])
    prod = [[sum(U[i][k] * B[k][j] for k in range(n)) for j in range(m)] for i in range(n)]
    return prod == Bred and abs(det_int(U)) == 1


def solve_row(B: Matrix, v: Sequence) -> Optional[List[Fraction]]:
    """
    Resuelve  x * B = v  con aritmetica exacta (Gauss-Jordan sobre Q).
    Devuelve x (lista de Fraction) o None si v no esta en el espacio generado por B.
    """
    n, m = len(B), len(B[0])
    A = [[Fraction(B[i][j]) for i in range(n)] + [Fraction(v[j])] for j in range(m)]
    r, piv_cols = 0, []
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None:
            raise ValueError("la base no tiene rango completo")
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [a * inv for a in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv_cols.append(c)
        r += 1
    if any(A[i][n] != 0 for i in range(r, m)):
        return None
    x = [Fraction(0)] * n
    for k, c in enumerate(piv_cols):
        x[c] = A[k][n]
    return x


def in_lattice(B: Matrix, v: Sequence) -> bool:
    """v pertenece a L(B) <=> sus coordenadas respecto a B son enteras."""
    x = solve_row(B, v)
    return x is not None and all(c.denominator == 1 for c in x)


def same_lattice_det(B: Matrix, Bred: Matrix) -> bool:
    """
    Comprueba misma celosia sin necesitar la matriz U (sirve para fpylll aunque yo no lo uso):
    si cada fila de Bred esta en L(B) y |det Bred| = |det B|, entonces L(Bred) es un
    sublattice de L(B) de indice 1, es decir, la misma celosia ([Nguyen, Lema 3]).
    """
    return all(in_lattice(B, r) for r in Bred) and abs(det_int(Bred)) == abs(det_int(B))


# Algoritmo LLL [LLL82] fig 1.13 de [smeets]
class Trace:
    """Registro de eventos para la memoria"""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled
        self.events: List[str] = []
        self.swaps = 0
        self.size_reductions = 0

    def log(self, msg: str):
        if self.enabled:
            self.events.append(msg)


def lll_reduce(
    B: Matrix,
    delta: Fraction = Fraction(3, 4),
    trace: Optional[Trace] = None,
) -> Tuple[Matrix, Matrix]:
    """
    Reduce la base B (filas) y devuelve (Bred, U) con Bred = U * B.

    Es la version iterativa de [LLL82]/[Smeets, fig. 1.13]:
      - mantiene mu_ij y B_i = ||b*_i||^2 y los actualiza tras cada intercambio
        (no recalcula la GSO entera);
      - k = indice del vector que se esta procesando (1-based en el articulo);
      - paso (1): size-reduce b_k respecto a b_{k-1} y comprueba Lovasz;
        si falla -> paso (2): intercambia b_{k-1} y b_k y retrocede (k-1);
        si se cumple -> size-reduce respecto al resto y avanza (k+1).
    delta en (1/4, 1]; el articulo usa 3/4. Para delta=1 no se conoce garantia
    polinomial [Nguyen], asi que se admite pero se avisa en la documentacion.
    """
    if not (Fraction(1, 4) < delta <= 1):
        raise ValueError("delta debe estar en (1/4, 1]")
    n = len(B)
    b = [row[:] for row in B]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    tr = trace or Trace(False)
    mu, Bn = gso(b)  # mu[i][j], Bn[i] = ||b*_i||^2 (exactos)

    def size_reduce(k: int, l: int):
        """Paso (*) del articulo: si |mu_kl| > 1/2, b_k -= r*b_l con r=redondeo(mu_kl)."""
        if abs(mu[k][l]) > Fraction(1, 2):
            # entero mas cercano; empates (|mu|=1/2 exacto) no llegan aqui
            r = (mu[k][l] + Fraction(1, 2)).__floor__()
            b[k] = [x - r * y for x, y in zip(b[k], b[l])]
            U[k] = [x - r * y for x, y in zip(U[k], U[l])]
            for j in range(l):
                mu[k][j] -= r * mu[l][j]
            mu[k][l] -= r
            tr.size_reductions += 1
            tr.log(f"  size-reduce: b{k+1} -= ({r})*b{l+1}")

    k = 1  # 0-based; corresponde a k=2 del articulo
    while k < n:
        size_reduce(k, k - 1)                                   # paso (1), l = k-1
        # condicion de Lovasz  B_k >= (delta - mu_{k,k-1}^2) B_{k-1}
        if Bn[k] < (delta - mu[k][k - 1] ** 2) * Bn[k - 1]:
            # paso (2): intercambio de b_{k-1} y b_k
            tr.swaps += 1
            tr.log(f"swap b{k} <-> b{k+1}   (Lovasz falla: B_k={float(Bn[k]):.4g} < "
                   f"{float((delta - mu[k][k-1]**2) * Bn[k-1]):.4g})")
            m = mu[k][k - 1]
            Bnew = Bn[k] + m * m * Bn[k - 1]
            mu[k][k - 1] = m * Bn[k - 1] / Bnew
            Bn[k] = Bn[k - 1] * Bn[k] / Bnew
            Bn[k - 1] = Bnew
            b[k - 1], b[k] = b[k], b[k - 1]
            U[k - 1], U[k] = U[k], U[k - 1]
            for j in range(k - 1):
                mu[k - 1][j], mu[k][j] = mu[k][j], mu[k - 1][j]
            # actualizacion de mu_{i,k-1}, mu_{i,k} para i > k  (formulas de [Cohen, 2.6.3])
            for i in range(k + 1, n):
                t = mu[i][k]
                mu[i][k] = mu[i][k - 1] - m * t
                mu[i][k - 1] = t + mu[k][k - 1] * mu[i][k]
            if k > 1:
                k -= 1
        else:
            for l in range(k - 2, -1, -1):                      # l = k-2, ..., 1 (1-based)
                size_reduce(k, l)
            k += 1
    return b, U


def shortest_basis_vector(B: Matrix) -> List[int]:
    """Vector de menor norma de una base (candidato a SVP tras LLL; no certifica el minimo)."""
    return min(B, key=norm2)
