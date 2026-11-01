"""Pruebas de rlwe.py: producto negaciclico, ataque en instancias generadas con el esquema LPR, ficheros reales."""
import os, random
from rlwe import (decrypt_bits, negacyclic, parse_rlwe_instance, polymul, solve_rlwe_instance, solve_rlwe_secret)

rng = random.Random(7)
q_of = {8: 257, 16: 1399}

# 1) el producto escolar coincide con la matriz negaciclica; x^n = -1
n, q = 8, 257
for _ in range(50):
    a = [rng.randrange(q) for _ in range(n)]; s = [rng.randrange(-2, 3) for _ in range(n)]
    A = negacyclic(a, n)
    assert polymul(a, s, q, n) == [sum(A[i][j] * s[j] for j in range(n)) % q for i in range(n)]
x = [0, 1] + [0] * (n - 2); xn1 = [0] * (n - 1) + [1]
assert polymul(x, xn1, q, n)[0] == q - 1                # x * x^(n-1) = x^n = -1
print("OK producto negaciclico")

# 2) instancias aleatorias del esquema (s,e en {-1,0,1}/{-2..2}); se recupera (s,e) y la palabra cifrada
fails = 0
for n in (8, 16):
    q = q_of[n]
    for t in range(6):
        a = [rng.randrange(0, 16) for _ in range(n)]
        s = [rng.randrange(-2, 3) for _ in range(n)]; e = [rng.randrange(-1, 2) for _ in range(n)]
        b = [(x + y) % q for x, y in zip(polymul(a, s, q, n), e)]
        r = solve_rlwe_secret(a, b, q)
        assert r["ok"] and r["consistent"], (n, t, r["note"])
        if r["secret"] != s:
            fails += 1                                   # solucion distinta pero valida y mas corta: permitido, se cuenta
        rr = [rng.randrange(-1, 2) for _ in range(n)]; e1 = [rng.randrange(-1, 2) for _ in range(n)]
        e2 = [rng.randrange(-1, 2) for _ in range(n)]; bits = [rng.randrange(2) for _ in range(n)]
        c1 = [(x + y) % q for x, y in zip(polymul(a, rr, q, n), e1)]
        c2 = [(x + y + bt * (q // 2)) % q for x, y, bt in zip(polymul(b, rr, q, n), e2, bits)]
        assert decrypt_bits(c1, c2, r["secret"], q)["bits"] == bits
print(f"OK 12 instancias aleatorias (secreto distinto del real en {fails}); descifrado correcto en todas")

# 3) ficheros reales n=16 y n=32 (si existen)
for f in ("rlwe_n16.txt", "rlwe_n32.txt"):
    p = os.path.join("instances", f)
    if os.path.exists(p):
        res = solve_rlwe_instance(parse_rlwe_instance(p), f)
        assert res["status"] == "OK_CERTIFIED" and res["detail"]["e_max"] <= 2 and len(res["bit"]) == res["k"], res
        print("OK", f, res["bit"])
