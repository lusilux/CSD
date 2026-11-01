"""
Pruebas del modulo LWE (lwe_io, lwe_attack, lwe_solve). Oraculos independientes:
  - el secreto/errores REALES con los que se generan las instancias aleatorias (modelo de Regev),
  - fuerza bruta sobre todos los s en Z_q^k (k pequeno),
  - los ejemplos 4 y 5 del PDF.
"""
import itertools
import os
import random
import tempfile
from fractions import Fraction

from lwe_attack import centered, decrypt_bit, pick_independent_rows, qary_lattice_basis, solve_secret
from lwe_io import is_lwe_file, parse_lwe_instance
from lwe_solve import solve_lwe_instance
from lll import det_int, in_lattice

# ---------------------------------------------------------------- 1) ejemplos 4 y 5 del PDF (q = 71)
S, E, q = [4, 70, 0, 3], [-3, 3, -2, 0, 2], 71
A = [[70, 6, 25, 15], [19, 51, 27, 32], [53, 55, 51, 6], [59, 46, 26, 51], [23, 5, 36, 45]]
b_pub = [32, 53, 31, 59, 11]
assert [(sum(a * x for a, x in zip(row, S)) + e) % q for row, e in zip(A, E)] == b_pub   # b = A s + e
u = [(A[0][j] + A[2][j]) % q for j in range(4)]                                          # Ej. 5: ecuaciones 1 y 3
assert u == [52, 61, 5, 21] and (b_pub[0] + b_pub[2]) % q == 63
c0 = decrypt_bit(u, 63, S, q)
assert c0["bit"] == 0 and c0["d"] == -5 == E[0] + E[2]                                    # el error acumulado es -5
c1 = decrypt_bit(u, (63 + q // 2) % q, S, q)                                             # 63+35 mod 71 = 27 (el PDF escribe 43)
assert (63 + q // 2) % q == 27 and c1["bit"] == 1
print("PDF Ej.4/5: b = A s + e, error acumulado d=-5 -> bit 0 y 63+35 mod 71 = 27 -> bit 1 "
      "(el PDF escribe 43: errata aritmetica)")

# ---------------------------------------------------------------- generador de Regev
random.seed(11)


def gen(k, m, q, B):
    s = [random.randrange(q) for _ in range(k)]
    A = [[random.randrange(q) for _ in range(k)] for _ in range(m)]
    e = [random.randint(-B, B) for _ in range(m)]
    b = [(sum(a * x for a, x in zip(row, s)) + ei) % q for row, ei in zip(A, e)]
    return s, A, e, b


def encrypt(A, b, q, bit):
    sub = [i for i in range(len(A)) if random.random() < 0.5] or [0]
    u = [sum(A[i][j] for i in sub) % q for j in range(len(A[0]))]
    v = (sum(b[i] for i in sub) + bit * (q // 2)) % q
    return u, v, sub


# ---------------------------------------------------------------- 2) instancias aleatorias: se recupera el secreto REAL
n_ok = n_noisy = 0
for k in range(2, 9):
    for qq in (71, 541, 1223):
        for B in (1, 5, 15):
            m = 2 * k + 8
            s, Am, e, b = gen(k, m, qq, B)
            r = solve_secret(Am, b, qq)
            if "secret" not in r:                # L tiene un vector mas corto que el error: solo admisible con ruido enorme
                assert (2 * B + 1) / qq > 0.1, (k, qq, B, r["note"])
                n_noisy += 1
                continue
            if r["spurious_log10"] > -3:         # solo se exige unicidad cuando es (casi) seguro
                continue
            assert r["ok"] and r["certified"] and r["cvp_agrees"], (k, qq, B)
            assert r["secret"] == s and r["errors_centered"] == e, (k, qq, B)
            for bit in (0, 1):
                u, v, sub = encrypt(Am, b, qq, bit)
                got = decrypt_bit(u, v, r["secret"], qq)["bit"]
                assert got == decrypt_bit(u, v, s, qq)["bit"]            # mismo resultado que con el secreto real
                if abs(sum(e[i] for i in sub)) < qq // 4:                # el esquema solo descifra bien si el error acumulado < q/4
                    assert got == bit, (k, qq, B)
            n_ok += 1
print(f"aleatorias OK: {n_ok} instancias, secreto y errores recuperados == los reales, bits 0 y 1 descifrados")

# ---------------------------------------------------------------- 3) fuerza bruta (oraculo independiente) con k = 2, 3
for k, qq in ((2, 71), (3, 71)):
    for m in (k + 1, k + 3, 2 * k + 2):
        s, Am, e, b = gen(k, m, qq, 3)
        best = min(sum(centered(b[i] - sum(a * x for a, x in zip(Am[i], t)), qq) ** 2 for i in range(m))
                   for t in itertools.product(range(qq), repeat=k))
        r = solve_secret(Am, b, qq)
        assert r["certified"] and r["err_norm2"] == best, (k, m, r["err_norm2"], best)
print("fuerza bruta en Z_q^k (k=2,3): la norma minima del error coincide con la del ataque")

# ---------------------------------------------------------------- 4) base q-aria y filas independientes
Am = [[3, 5, 7], [1, 1, 1], [2, 4, 6], [0, 1, 9], [8, 2, 2]]
L, perm = qary_lattice_basis(Am, 71)
assert abs(det_int(L)) == 71 ** (5 - 3)
for i in range(5):                                       # A s mod q esta en la celosia (en el orden permutado)
    pass
sv = [10, 20, 30]
y = [sum(a * x for a, x in zip(Am[i], sv)) % 71 for i in perm]
assert in_lattice(L, y)
assert len(pick_independent_rows(Am, 71)) == 3
try:
    pick_independent_rows([[1, 2], [2, 4], [3, 6]], 71)
    raise AssertionError("deberia fallar con rango < k")
except ValueError:
    pass
print("celosia q-aria: det = q^(m-k), contiene a A s mod q; deteccion de rango deficiente")

# ---------------------------------------------------------------- 5) ficheros reales de la asignatura
for name, bit in (("instances/lwe_t3_k7_m18.txt", 0), ("instances/lwe_t3_k10_m30.txt", 1)):
    if not os.path.exists(name):
        continue
    inst = parse_lwe_instance(name)
    assert is_lwe_file(name)
    res = solve_lwe_instance(inst, os.path.basename(name))
    assert res["status"] == "OK_CERTIFIED" and res["bit"] == bit, (name, res["status"], res["bit"])
    d = res["detail"]
    assert d["err_max"] <= 15 and d["spurious_log10"] < -5
    assert d["consistent"] and d["cvp_agrees"]
print("ficheros LWE de la asignatura: certificados, consistentes, con CVP directo coincidente")

# ---------------------------------------------------------------- 6) menos muestras, presupuesto, errores de formato
inst = parse_lwe_instance("instances/lwe_t3_k10_m30.txt")
r24 = solve_secret(inst["A"], inst["b"], inst["q"], samples=24)
assert r24["ok"] and r24["secret"] == solve_secret(inst["A"], inst["b"], inst["q"])["secret"]
rb = solve_lwe_instance(inst, "budget", max_nodes=1)
assert rb["status"] in ("UNCERTIFIED_BUDGET", "FAILED_NO_SECRET", "FAILED_VERIFICATION") and not rb["certified"]
with tempfile.TemporaryDirectory() as td:
    p = os.path.join(td, "mal.txt")
    open(p, "w").write("rotura\nN/A\n3\n71\n{1,2,3},\n{4,5},\n{1,2}\nN/D\n")
    try:
        parse_lwe_instance(p)
        raise AssertionError("deberia rechazar filas con distinto numero de coeficientes")
    except ValueError:
        pass
    open(p, "w").write("SVP\n0 0\n2\n1 0\n0 1\nN/A\nN/D\n")
    assert not is_lwe_file(p)
print("muestras reducidas (24 de 30), presupuesto agotado y ficheros mal formados: OK")
print("todas las pruebas LWE OK")
