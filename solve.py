"""
solve.py 
Uso:  python solve.py instances/*.txt [--delta 0.99] [--backend auto|exact|fpylll]
                                      [--max-nodes N] [--seconds S] [--csv salida.csv] [--trace]
"""
import argparse
import csv
import datetime
import glob
import json
import os
import platform
import statistics
import sys
import time
from fractions import Fraction
from typing import Dict, List, Optional

from babai import babai_nearest_plane, dist2
from enumeration import solve_cvp, solve_svp
from instance_io import parse_lattice_instance
from lll import (Matrix, Trace, det_int, in_lattice, is_lll_reduced, lll_reduce, norm2,
                 same_lattice, same_lattice_det)
from verify import coefficients, gaussian_heuristic2, verify_cvp, verify_svp
from lwe_io import is_lwe_file, parse_lwe_instance
from rlwe import is_rlwe_file

EXACT_MAX_DIM = 30   # por encima, backend='auto' usa fpylll. Esto no se usa en los ejemplos de clase


def reduce_basis(B: Matrix, delta: Fraction, backend: str = "auto", trace: Optional[Trace] = None,
                 stats: Optional[dict] = None):
    """Devuelve (Bred, backend_usado). Comprueba siempre que Bred genera la misma celosia.
    Si se pasa `stats` se rellena con t_lll, t_check y swaps (swaps: solo backend exacto)."""
    n = len(B)
    if backend == "auto":
        backend = "exact" if n <= EXACT_MAX_DIM else "fpylll"
    if backend == "exact":
        tr = trace or Trace(False)                       # Trace(False) solo cuenta intercambios
        t = time.perf_counter()
        Bred, U = lll_reduce(B, delta, tr)
        t_lll = time.perf_counter() - t
        t = time.perf_counter()
        assert is_lll_reduced(Bred, delta), "salida no LLL-reducida"
        assert same_lattice(B, U, Bred), "la base reducida no genera la misma celosia"
        if stats is not None:
            stats.update(t_lll=t_lll, t_check=time.perf_counter() - t, swaps=tr.swaps)
        return Bred, "exact"
    if backend == "fpylll":
        from fpylll import LLL, IntegerMatrix   # Ahora mismo esto no funciona
        t = time.perf_counter()
        A = IntegerMatrix.from_matrix(B)
        LLL.reduction(A, delta=float(delta), eta=0.501)
        Bred = [[A[i, j] for j in range(A.ncols)] for i in range(A.nrows)]
        t_lll = time.perf_counter() - t
        t = time.perf_counter()
        assert same_lattice_det(B, Bred), "la base reducida no genera la misma celosia"
        if stats is not None:
            stats.update(t_lll=t_lll, t_check=time.perf_counter() - t, swaps=None)   # fpylll no expone swaps
        return Bred, "fpylll"
    raise ValueError(f"backend desconocido: {backend}")


def _canonical_sign(v: List[int]) -> List[int]:
    """Entre v y -v se devuelve el que tiene la primera coordenada no nula positiva."""
    for a in v:
        if a != 0:
            return list(v) if a > 0 else [-b for b in v]
    return list(v)

def celosia_modular(B, n):
    """Base de L(B) + n*Z^k: filas de B mas los vectores n*e_j, eliminacion entera por columnas."""
    k = len(B[0])
    pool = [[x % n for x in row] for row in B]
    pivotes = []
    for j in range(k):
        pool.append([n if c == j else 0 for c in range(k)])   # entra n*e_j
        while True:                                            # Euclides en la columna j
            nz = [r for r in pool if r[j] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: abs(r[j]))
            p = nz[0]
            for r in nz[1:]:
                q = r[j] // p[j]
                for c in range(k):
                    r[c] -= q * p[c]
        piv = next(r for r in pool if r[j] != 0)
        pool.remove(piv)
        if piv[j] < 0:
            piv = [-x for x in piv]
        pivotes.append(piv)
        for r in pool:                                         # valido: n*e_c sigue en la celosia
            for c in range(j + 1, k):
                r[c] %= n
    return pivotes

def solve_instance(inst: Dict, name: str = "", delta: Fraction = Fraction(99, 100),
                   backend: str = "auto", max_nodes=None, max_seconds=None,
                   trace: Optional[Trace] = None) -> Dict:
    B, n, problem = inst["basis"], inst["dim"], inst["problem"]
    out = dict(name=name, problem=problem, dim=n, status="", answer=None, value=None,
               coeffs=None, certified=False, verified=False, backend="", t_reduce=0.0,
               t_solve=0.0, nodes=0, source="", initial_value=None, note="",
               value_kind="norma^2" if problem == "SVP" else "distancia^2",
               modulus="N/A", reading="N/A")
    if inst["modulus"] is not None:
        B = celosia_modular(B, inst["modulus"])
        out["note"] = f"celosia L(B)+{inst['modulus']}*Z^{n} (coef. respecto a la nueva base)"
        out["modulus"], out["reading"] = inst["modulus"], f"L(B)+nZ^k, n={inst['modulus']}"
    if problem == "CVP" and len(inst["target"]) != n:
        out.update(status="ERROR", note="la dimension del objetivo no coincide con la de la base")
        return out

    if det_int(B) == 0:
        out.update(status="ERROR", note="la base no tiene rango completo (det = 0)")
        return out

    t0 = time.perf_counter()
    st: dict = {}
    Bred, used = reduce_basis(B, delta, backend, trace, st)
    out["t_reduce"], out["backend"] = time.perf_counter() - t0, used
    out.update(t_lll=st["t_lll"], t_check=st["t_check"], swaps=st["swaps"],
               lll_first2=min(norm2(r) for r in Bred), gh2=gaussian_heuristic2(B))

    if problem == "SVP":
        r = solve_svp(Bred, max_nodes, max_seconds)
        v = _canonical_sign(r.vector)
        check = verify_svp(B, v)
        out["note"] = out["note"] + " | " + f"gauss_heuristic^2~{gaussian_heuristic2(B):.1f}"
    elif problem == "CVP":
        t = inst["target"]
        r = solve_cvp(Bred, t, max_nodes, max_seconds)
        v = r.vector
        check = verify_cvp(B, t, v)
    else:
        out.update(status="ERROR", note=f"problema desconocido: {problem}")
        return out

    out.update(answer=v, value=r.dist2, coeffs=coefficients(B, v), certified=r.certified,
               verified=check["ok"], t_solve=r.seconds, nodes=r.nodes, source=r.source,
               initial_value=r.initial_dist2)
    out["status"] = ("OK_CERTIFIED" if r.certified else "UNCERTIFIED_BUDGET") if check["ok"] else "FAILED_VERIFICATION"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--delta", default="0.99")
    ap.add_argument("--backend", default="auto", choices=["auto", "exact", "fpylll"])
    ap.add_argument("--max-nodes", type=int, default=None)
    ap.add_argument("--seconds", type=float, default=None)
    ap.add_argument("--samples", type=int, default=None, help="LWE: usar solo las m' primeras ecuaciones")
    ap.add_argument("--embed-M", type=int, default=1, help="LWE: factor M de la incrustacion de Kannan")
    ap.add_argument("--no-crosscheck", action="store_true", help="LWE: omitir la verificacion CVP directa")
    ap.add_argument("--repeat", type=int, default=1,
                    help="repite las instancias rapidas (t_total < 0.1 s) N veces y guarda la mediana de los tiempos")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--trace", action="store_true", help="imprime la traza de LLL")
    a = ap.parse_args(argv)
    delta = Fraction(a.delta)

    # PowerShell/cmd no expanden comodines (instances\\*.txt): lo hacemos aqui
    paths = []
    for p in a.files:
        paths.extend(sorted(glob.glob(p)) or [p])

    rows = []

    def timed(fn):
        """Ejecuta fn() midiendo t_total (todo el trabajo, con comprobaciones). Si es muy rapida y --repeat > 1,
        repite y deja la mediana de t_lll, t_check, t_solve y t_total (el resultado es el de la ultima)."""
        def once():
            t = time.perf_counter()
            r = fn()
            r["t_total"] = time.perf_counter() - t
            return r
        res = once()
        runs = [res]
        if a.repeat > 1 and res["t_total"] < 0.1:
            runs += [once() for _ in range(a.repeat - 1)]
            res = runs[-1]
            for key in ("t_lll", "t_check", "t_solve", "t_total"):
                vals = [r[key] for r in runs if r.get(key) is not None]
                if vals:
                    res[key] = statistics.median(vals)
        res["repeats"] = len(runs)
        res["t_enum"] = res.get("t_solve")
        res["date"] = datetime.datetime.now().isoformat(timespec="seconds")
        if res.get("answer") is not None:
            res["answer_text"] = json.dumps(res["answer"])
        return res

    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        if is_rlwe_file(path):                     # ficheros Ring-LWE (anillo x^n+1)
            from rlwe import parse_rlwe_instance, print_rlwe_result, solve_rlwe_instance
            res = timed(lambda: solve_rlwe_instance(parse_rlwe_instance(path), name, delta, a.backend,
                                                    a.max_nodes, a.seconds, a.embed_M))
            rows.append(res)
            print_rlwe_result(res)
            continue
        if is_lwe_file(path):                      # ficheros LWE (rotura / criptoanalisis)
            from lwe_solve import print_lwe_result, solve_lwe_instance
            res = timed(lambda: solve_lwe_instance(parse_lwe_instance(path), name, delta, a.backend, a.max_nodes,
                                                   a.seconds, a.samples, a.embed_M, not a.no_crosscheck))
            rows.append(res)
            print_lwe_result(res)
            continue
        tr = Trace(True) if a.trace else None
        res = timed(lambda: solve_instance(parse_lattice_instance(path), name, delta, a.backend,
                                           a.max_nodes, a.seconds, tr))
        rows.append(res)
        print(f"\n== {name}  {res['problem']} dim={res['dim']}  -> {res['status']}")
        if res["status"] in ("OK_CERTIFIED", "UNCERTIFIED_BUDGET", "FAILED_VERIFICATION"):
            kind = "norma^2" if res["problem"] == "SVP" else "distancia^2"
            print(f"   respuesta : {res['answer']}   ({kind} = {res['value']}; "
                  f"partida = {res['initial_value']}, origen: {res['source']})")
            print(f"   coef. base original: {res['coeffs']}")
            print(f"   backend={res['backend']}  t_LLL={res['t_reduce']:.4f}s  "
                  f"t_enum={res['t_solve']:.4f}s  nodos={res['nodes']}  "
                  f"OPTIMO_CERTIFICADO={res['certified']}  comprobacion_basica(pertenece a la celosia{' + Hermite' if res['problem']=='SVP' else ''})={res['verified']}")
            if res["note"]:
                print("   nota:", res["note"])
        else:
            print("   ", res["note"])
        if tr:
            print("   traza LLL:\n    " + "\n    ".join(tr.events[:40]) +
                  ("\n    ..." if len(tr.events) > 40 else ""))

    if a.csv:
        cols = ["name", "problem", "dim", "status", "value", "value_kind", "certified", "verified", "backend",
                "t_lll", "t_check", "t_enum", "t_total", "repeats", "swaps", "nodes", "source",
                "modulus", "reading", "lll_first2", "gh2",                       # celosias
                "k", "m", "q", "bit", "s_max", "e_max", "gap", "spurious_log10",  # LWE / Ring-LWE
                "answer_text", "date",
                "t_reduce", "t_solve"]                                           # columnas antiguas (compatibilidad)
        # Se ACUMULA entre ejecuciones: si el CSV ya existe se conservan sus filas y se
        # reemplazan solo las instancias con el mismo nombre (sin duplicados al repetir).
        old = {}
        if os.path.exists(a.csv):
            with open(a.csv, newline="", encoding="utf-8") as f:
                old = {r["name"]: r for r in csv.DictReader(f)}
        for r in rows:
            old[r["name"]] = {k: r.get(k) for k in cols}
        with open(a.csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(old.values())
        print(f"\nCSV actualizado: {a.csv} ({len(old)} instancias)")
        meta = dict(fecha=datetime.datetime.now().isoformat(timespec="seconds"), sistema=platform.platform(),
                    cpu=platform.processor() or platform.machine(), nucleos=os.cpu_count(),
                    python=platform.python_version(), delta=a.delta, backend_pedido=a.backend, repeat=a.repeat)
        try:
            import fpylll
            meta["fpylll"] = fpylll.__version__
        except ImportError:
            meta["fpylll"] = None
        with open(os.path.splitext(a.csv)[0] + "_equipo.json", "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2, ensure_ascii=False)
    return rows


if __name__ == "__main__":
    sys.exit(0 if main() is not None else 1)
