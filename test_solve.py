"""
Pruebas del pipeline SVP/CVP: Babai, enumeracion y solve.py.
Oraculos independientes: fpylll (enumeracion/CVP de fplll) y fuerza bruta en caja (dim pequena).
"""
import glob
import os
import random
from fractions import Fraction

from fpylll import CVP, GSO, LLL, Enumeration, IntegerMatrix

from babai import babai_nearest_plane, babai_round, dist2
from enumeration import _zigzag, enumerate_vectors, solve_cvp, solve_svp
from instance_io import parse_lattice_instance
from lwe_io import is_lwe_file
from rlwe import is_rlwe_file
from lll import det_int, in_lattice, lll_reduce, norm2, same_lattice_det, solve_row
from solve import solve_instance
from verify import brute_force_cvp, brute_force_svp, coefficients


# ---------------------------------------------------------------- oraculos (fpylll)
def oracle_svp_norm2(B):
    n = len(B)
    A = IntegerMatrix.from_matrix(B)
    M = GSO.Mat(A)
    M.update_gso()
    LLL.Reduction(M)()
    _, coeffs = Enumeration(M, nr_solutions=1).enumerate(0, n, M.get_r(0, 0) + 0.5, 0)[0]
    v = [int(round(sum(c * A[i, j] for i, c in enumerate(coeffs)))) for j in range(n)]
    return norm2(v)


def oracle_cvp_dist2(B, t):
    A = IntegerMatrix.from_matrix(B)
    LLL.reduction(A)          # fplll CVP falla con algunas bases sin reducir (ReductionError 5)
    p = list(CVP.closest_vector(A, list(t)))
    return dist2(p, t)


# ---------------------------------------------------------------- 1) ejemplos del PDF (GGH)
B_priv, B_pub, U = [[3, 1], [-1, 4]], [[5, 19], [3, 14]], [[3, 4], [2, 3]]
Uinv = [[3, -4], [-2, 3]]
c = [43, 176]                                    # Ejemplo 2: x=(7,3), e=(-1,1)
v, y = babai_round(B_priv, c)
assert y == [27, 37], y                          # Ejemplo 3: y = Round(c B^-1) = (27, 37)
x = [sum(y[i] * Uinv[i][j] for i in range(2)) for j in range(2)]
assert x == [7, 3], x                            # x = y U^-1 = (7, 3)
c2 = [42, 177]                                   # e' = (-2, 2): el PDF dice que falla
v2, y2 = babai_round(B_priv, c2)
x2 = [sum(y2[i] * Uinv[i][j] for i in range(2)) for j in range(2)]
assert x2 != [7, 3], "Babai round deberia fallar con e'=(-2,2)"
r = solve_cvp(B_pub, c2)                         # CVP EXACTO sobre la base publica
assert r.certified and r.dist2 == oracle_cvp_dist2(B_pub, c2) == 5
assert dist2([44, 175], c2) == 8 and r.vector == [43, 179]
print(f"PDF Ej.3: con e'=(-2,2) Babai da x={x2}; el punto MAS cercano a c' es {r.vector} "
      f"(dist2={r.dist2}) y no el original (44,175) (dist2=8): la ambiguedad es real, no un fallo de Babai.")

# ---------------------------------------------------------------- 1b) orden de _zigzag (la poda 'break' depende de el)
from enumeration import round_frac
for q in range(1, 13):
    for p in range(-60, 61):
        c = Fraction(p, q)
        x0 = round_frac(c)
        z = _zigzag(x0, 1 if c - x0 >= 0 else -1)
        ds = [abs(next(z) - c) for _ in range(12)]
        assert all(a <= b for a, b in zip(ds, ds[1:])), (c, ds)
print("zigzag: distancias no decrecientes para centros p/q (empates incluidos)")

# ---------------------------------------------------------------- 2) algebra exacta
assert solve_row(B_pub, [43, 179]) == [5, 6]
assert not in_lattice(B_pub, [43, 180])                      # punto fuera de la celosia
assert in_lattice(B_pub, [44, 175]) and coefficients(B_pub, [44, 175]) == [7, 3]

# ---------------------------------------------------------------- 3) todas las instancias del usuario
checked = 0
for path in sorted(glob.glob("instances/*.txt")):
    if is_lwe_file(path) or is_rlwe_file(path):   # los ficheros LWE / Ring-LWE los prueban test_lwe.py y test_rlwe.py
        continue
    name = os.path.splitext(os.path.basename(path))[0]
    inst = parse_lattice_instance(path)
    res = solve_instance(inst, name)
    if inst["modulus"] is not None:
        # lectura q-aria L(B)+n*Z^k (tu solve.py ya la aplica): con B invertible mod n sale Z^k
        # (SVP: norma^2 1; CVP: distancia 0). No hay oraculo fpylll para esta lectura.
        assert res["status"] == "OK_CERTIFIED", (name, res["status"])
        continue
    assert res["status"] == "OK_CERTIFIED", (name, res["status"])
    ref = (oracle_svp_norm2(inst["basis"]) if inst["problem"] == "SVP"
           else oracle_cvp_dist2(inst["basis"], inst["target"]))
    assert res["value"] == ref, (name, res["value"], ref)
    checked += 1
print(f"instancias: {checked} resueltas, certificadas, verificadas y iguales al oraculo fpylll")

# ---------------------------------------------------------------- 4) bases aleatorias (tipo GGH, B' = U G)
random.seed(7)


def random_unimodular(n, steps=60, mag=3):
    M = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(steps):
        i, j = random.sample(range(n), 2)
        k = random.randint(-mag, mag)
        M[i] = [a + k * b for a, b in zip(M[i], M[j])]
    return M


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


improved_over_shortest_row = babai_wrong = far_cases = 0
for n in range(2, 10):
    for _ in range(6):
        while True:
            G = [[random.randint(-9, 9) for _ in range(n)] for _ in range(n)]
            if det_int(G) != 0:
                break
        B = matmul(random_unimodular(n), G)
        for delta in (Fraction(3, 4), Fraction(99, 100)):
            Bred, _ = lll_reduce(B, delta)
            # SVP
            s = solve_svp(Bred)
            assert s.certified and in_lattice(B, s.vector) and any(s.vector)
            assert s.dist2 == oracle_svp_norm2(B), (n, s.dist2)
            if n <= 4:
                assert s.dist2 == brute_force_svp(Bred, 6)[0]
            improved_over_shortest_row += s.source == "enumeracion"
            # CVP: objetivo cercano y objetivo lejano
            for scale in (3, 400):
                t = [random.randint(-scale, scale) for _ in range(n)]
                q = solve_cvp(Bred, t)
                assert q.certified and in_lattice(B, q.vector)
                assert q.dist2 == oracle_cvp_dist2(B, t), (n, scale)
                if n <= 4:
                    assert q.dist2 == brute_force_cvp(Bred, t, 4)[0]
                far_cases += scale == 400
                babai_wrong += q.source == "enumeracion"
print(f"aleatorias OK. La enumeracion mejoro a la fila mas corta en {improved_over_shortest_row} SVP "
      f"y a Babai en {babai_wrong} CVP (de {far_cases * 2} CVP; la mitad con objetivo lejano)")

# ---------------------------------------------------------------- 5) backend fpylll
for n in (5, 12, 20):
    G = [[random.randint(-9, 9) for _ in range(n)] for _ in range(n)]
    while det_int(G) == 0:
        G = [[random.randint(-9, 9) for _ in range(n)] for _ in range(n)]
    B = matmul(random_unimodular(n), G)
    inst = dict(problem="SVP", target=[0] * n, dim=n, basis=B, modulus=None, solution="N/D")
    a = solve_instance(inst, f"fp{n}", backend="fpylll")
    b = solve_instance(inst, f"ex{n}", backend="exact")
    assert a["status"] == b["status"] == "OK_CERTIFIED" and a["value"] == b["value"], (n, a["value"], b["value"])
print("backend fpylll == backend exacto (SVP, dims 5/12/20)")

# ---------------------------------------------------------------- 6) presupuesto agotado
inst = parse_lattice_instance("instances/svp_dim15.txt")
Bred, _ = lll_reduce(inst["basis"], Fraction(99, 100))
r0 = enumerate_vectors(Bred, None, radius2=10**9, max_nodes=3)   # radio enorme -> arbol grande
assert not r0.certified and r0.nodes > 3 - 1
res = solve_instance(dict(inst, basis=inst["basis"]), "budget", max_nodes=1)
assert res["status"] == "UNCERTIFIED_BUDGET" and res["verified"] and not res["certified"], res["status"]
print("presupuesto agotado -> se informa 'no certificado', sin excepcion y con respuesta valida")

# ---------------------------------------------------------------- 7) nunca devuelve el vector nulo
for B in ([[1, 0], [0, 1]], [[2, 0, 0], [0, 3, 0], [0, 0, 5]]):
    s = solve_svp(B)
    assert any(s.vector) and s.certified
print("todas las pruebas OK")
