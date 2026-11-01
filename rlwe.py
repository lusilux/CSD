"""
rlwe.py -- Ring-LWE (anillo R_q = Z_q[x]/(x^n+1)) atacado con los MISMOS bloques que LWE/SVP:
LLL + enumeracion sobre una celosia, via incrustacion de Kannan.

Modelo (enunciado): polinomios a, b de grado < n con  b = a*s + e  en R_q  (s, e pequenos).
  cifrado de n bits (Lyubashevsky-Peikert-Regev 2010): c1 = a*r + e1,  c2 = b*r + e2 + bits*floor(q/2).
  descifrado: d = c2 - c1*s = e2 + e*r - e1*s + bits*floor(q/2)  (en R_q); cada coeficiente
  |d_i| pequeno -> bit 0 ; |d_i - q/2| pequeno -> bit 1.   (LPR10 de memoria; VERIFICAR antes de citar)

Reduccion a celosia (analoga a LWE, ver lwe_attack.py):
  El producto a*s (mod x^n+1) es lineal en s: matriz NEGACICLICA A, A[i][j] = a[i-j] si i>=j,
  -a[n+i-j] si i<j  (porque x^n = -1).  Entonces b = A s + e (mod q), un LWE con k=m=n.
  Con m = k las ecuaciones solas NO determinan s (A es invertible y cualquier b es alcanzable);
  lo que lo determina es que s Y e son pequenos.  Por eso se usa la celosia de dimension 2n
      L = { (s, A s mod q) }  base:  [ I | A^T ]  y  [ 0 | q I ]       (det = q^n)
  y se incrusta el objetivo (0, b): fila (0, b, M).  El vector (-s, e, M) es muy corto -> SVP unico.
  (Es la 'forma normal' del primal attack con (s,e) conjuntos.)  Se resuelve con lll.py + enumeration.py
  (o fpylll para LLL cuando la dimension es grande, igual que en lwe_attack.py).

Que se certifica: la enumeracion completa demuestra que (-s,e,M) es el vector no nulo mas corto de la
celosia incrustada (norma minima de (s,e)); que sea el par 'real' del cifrador es probable por el
hueco, no demostrado.  Comprobacion independiente: b == a*s + e en R_q con un producto polinomico
hecho aparte (no con la matriz).
"""
import re
import time
from fractions import Fraction
from typing import Dict, List, Sequence

from enumeration import solve_svp
import math
from lwe_attack import _reduce, centered
from lwe_io import _data_lines


def is_rlwe_file(path: str) -> bool:
    lines = _data_lines(path)
    return bool(lines) and "ring" in lines[0].lower()


def _ints(s):
    return [int(x) for x in re.findall(r"-?\d+", s)]


def parse_rlwe_instance(path: str) -> Dict:
    L = _data_lines(path)
    cipher = None
    if L[1].strip().upper() != "N/A":
        cipher = (_ints(L[1]), _ints(L[2]))
        L = L[:1] + [L[1]] + L[3:]                       # quita la 2a linea del criptograma
    n, q = int(L[2]), int(L[3])
    a, b = _ints(L[4]), _ints(L[5])
    if not (len(a) == len(b) == n) or (cipher and not (len(cipher[0]) == len(cipher[1]) == n)):
        raise ValueError(f"{path}: se esperaban polinomios con n={n} coeficientes")
    return dict(problem="RLWE", task="criptoanalisis" if cipher else "rotura", cipher=cipher, n=n, q=q,
                a=a, b=b, solution=L[6] if len(L) > 6 else "N/D")


def polymul(a: Sequence[int], s: Sequence[int], q: int, n: int) -> List[int]:
    """a*s en Z_q[x]/(x^n+1), producto escolar (independiente de la matriz negaciclica)."""
    r = [0] * n
    for i, ai in enumerate(a):
        for j, sj in enumerate(s):
            k = i + j
            if k < n:
                r[k] += ai * sj
            else:
                r[k - n] -= ai * sj
    return [x % q for x in r]


def negacyclic(a: Sequence[int], n: int) -> List[List[int]]:
    return [[a[i - j] if i >= j else -a[n + i - j] for j in range(n)] for i in range(n)]


def solve_rlwe_secret(a, b, q, M=1, delta=Fraction(99, 100), backend="auto", max_nodes=None, max_seconds=None):
    n = len(a)
    A = negacyclic(a, n)
    basis = [[int(i == j) for i in range(n)] + [A[r][j] % q for r in range(n)] for j in range(n)]
    basis += [[0] * n + [q if r == i else 0 for r in range(n)] for i in range(n)]
    E = [row + [0] for row in basis] + [[0] * n + [x % q for x in b] + [M]]
    t0 = time.perf_counter()
    st: dict = {}
    Ered, used = _reduce(E, delta, backend, st)
    out = dict(ok=False, n=n, q=q, dim=2 * n + 1, backend=used, t_reduce=time.perf_counter() - t0, note="",
               t_lll=st["t_lll"], t_check=st["t_check"], swaps=st["swaps"])
    r = solve_svp(Ered, max_nodes, max_seconds)
    w = list(r.vector)
    if w[-1] < 0:
        w = [-x for x in w]
    out.update(t_solve=r.seconds, nodes=r.nodes, certified=r.certified, source=r.source, value=sum(x * x for x in w))
    if w[-1] != M:
        out["note"] = "la incrustacion no aisla (s,e): el vector mas corto no lleva la coordenada M"
        return out
    s, e = [(-x) % q for x in w[:n]], [x % q for x in w[n:2 * n]]
    check = polymul(a, s, q, n)
    consistent = all((b[i] - check[i] - e[i]) % q == 0 for i in range(n))
    sc, ec = [centered(x, q) for x in s], [centered(x, q) for x in e]
    out.update(secret=sc, errors=ec, consistent=consistent, s_max=max(map(abs, sc)), e_max=max(map(abs, ec)),
               norm2=sum(x * x for x in sc + ec))
    # hueco (heuristica gaussiana de la celosia de dim 2n y det q^n) y soluciones espurias esperadas:
    # pares (s,e) con todas las entradas en [-B,B] que cumplen b = a*s+e:  (2B+1)^(2n) / q^n
    gh = math.sqrt(2 * n / (2 * math.pi * math.e)) * math.sqrt(q)
    B_ = max(out["s_max"], out["e_max"])
    out["gap"] = gh / math.sqrt(out["norm2"]) if out["norm2"] else float("inf")
    out["spurious_log10"] = 2 * n * math.log10(2 * B_ + 1) - n * math.log10(q)
    out["ok"] = consistent and out["certified"]
    return out


def decrypt_bits(c1, c2, s, q) -> Dict:
    n = len(s)
    c1s = polymul(c1, s, q, n)
    half, bits, margins = q // 2, [], []
    for i in range(n):
        d = c2[i] - c1s[i]
        d0, d1 = abs(centered(d, q)), abs(centered(d - half, q))
        bits.append(0 if d0 < d1 else 1)
        margins.append(min(d0, d1))
    return dict(bits=bits, max_noise=max(margins), threshold=q // 4)


def solve_rlwe_instance(inst: Dict, name="", delta=Fraction(99, 100), backend="auto", max_nodes=None,
                        max_seconds=None, M=1) -> Dict:
    n, q = inst["n"], inst["q"]
    row = dict(name=name, problem="RLWE-" + inst["task"], dim=2 * n + 1, status="", answer=None, value=None,
               certified=False, verified=False, backend="", t_reduce=0.0, t_solve=0.0, nodes=0, source="",
               k=n, m=n, q=q, bit=None, note="", detail={})
    try:
        r = solve_rlwe_secret(inst["a"], inst["b"], q, M, delta, backend, max_nodes, max_seconds)
    except ValueError as ex:
        row.update(status="ERROR", note=str(ex))
        return row
    row.update(backend=r["backend"], t_reduce=r["t_reduce"], t_solve=r.get("t_solve", 0.0), nodes=r.get("nodes", 0),
               source=r.get("source", ""), certified=r.get("certified", False), value=r.get("norm2"),
               detail=r, note=r["note"])
    row.update(t_lll=r.get("t_lll"), t_check=r.get("t_check"), swaps=r.get("swaps"), value_kind="||(s,e)||^2")
    if "secret" not in r:
        row["status"] = "FAILED_NO_SECRET"
        return row
    row.update(e_max=r["e_max"], s_max=r["s_max"], gap=r["gap"], spurious_log10=r["spurious_log10"])
    row["answer"] = r["secret"]
    row["verified"] = r["consistent"]
    if inst["cipher"] is not None:
        d = decrypt_bits(inst["cipher"][0], inst["cipher"][1], r["secret"], q)
        row["bit"] = "".join(map(str, d["bits"]))
        r["decrypt"] = d
        if d["max_noise"] >= d["threshold"]:
            row["note"] += " AVISO: ruido de descifrado >= q/4, la palabra puede ser incorrecta."
    row["status"] = ("OK_CERTIFIED" if r["certified"] else "UNCERTIFIED_BUDGET") if row["verified"] else "FAILED_VERIFICATION"
    return row


def print_rlwe_result(res: Dict) -> None:
    print(f"\n== {res['name']}  {res['problem']}  n={res['k']} q={res['q']}  -> {res['status']}")
    if res["answer"] is None:
        print("   ", res["note"])
        return
    d = res["detail"]
    print(f"   secreto s (coef. x^0..x^(n-1)) = {res['answer']}")
    print(f"   errores e = {d['errors']}   (max|s|={d['s_max']}, max|e|={d['e_max']}, ||(s,e)||^2={d['norm2']})")
    if res["bit"] is not None:
        c = d["decrypt"]
        print(f"   PALABRA DE BITS (bit i = coef. x^i) = {res['bit']}   (ruido max por coef. = {c['max_noise']}, "
              f"umbral q/4 = {c['threshold']})")
    print(f"   celosia incrustada dim={d['dim']}  backend={res['backend']}  t_LLL={res['t_reduce']:.2f}s  "
          f"t_enum={res['t_solve']:.2f}s  nodos={res['nodes']}  OPTIMO_CERTIFICADO={res['certified']}")
    print(f"   comprobacion: b = a*s + e en R_q (producto polinomico independiente) = {d['consistent']}")
    if res["note"]:
        print("   nota:", res["note"])
