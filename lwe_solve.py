"""
lwe_solve.py -- Interfaz de las tareas LWE (rotura / criptoanalisis) para solve.py.
La logica matematica esta en lwe_attack.py; aqui solo se orquesta y se formatea la salida.
"""
from fractions import Fraction
from typing import Dict, Optional

from lwe_attack import centered, decrypt_bit, solve_secret
from lwe_io import parse_lwe_instance


def solve_lwe_instance(inst: Dict, name: str = "", delta: Fraction = Fraction(99, 100), backend: str = "auto",
                       max_nodes=None, max_seconds=None, samples: Optional[int] = None, M: int = 1,
                       cross_check: bool = True) -> Dict:
    k, m, q = inst["k"], inst["m"], inst["q"]
    row = dict(name=name, problem="LWE-" + inst["task"], dim=(samples or m) + 1, status="", answer=None,
               value=None, certified=False, verified=False, backend="", t_reduce=0.0, t_solve=0.0,
               nodes=0, source="", k=k, m=m, q=q, bit=None, note="", detail={})
    try:
        r = solve_secret(inst["A"], inst["b"], q, M, delta, backend, samples, cross_check, max_nodes, max_seconds)
    except ValueError as ex:
        row.update(status="ERROR", note=str(ex))
        return row
    row.update(backend=r.get("backend", ""), t_reduce=r.get("t_reduce", 0.0), t_solve=r.get("t_solve", 0.0),
               nodes=r.get("nodes", 0), source=r.get("source", ""), certified=r.get("certified", False),
               value=r.get("err_norm2"), detail=r, note=r.get("note", ""))
    row.update(t_lll=r.get("t_lll"), t_check=r.get("t_check"), swaps=r.get("swaps"), value_kind="||e||^2")
    if "secret" not in r:
        row["status"] = "FAILED_NO_SECRET"
        return row
    row.update(e_max=r["err_max"], s_max=max(abs(centered(x, q)) for x in r["secret"]), gap=r["gap"],
               spurious_log10=r["spurious_log10"])
    row["answer"] = r["secret"]
    row["verified"] = bool(r["consistent"]) and r.get("cvp_agrees", True)
    if inst["cipher"] is not None:
        d = decrypt_bit(inst["cipher"][0], inst["cipher"][1], r["secret"], q)
        row["bit"], row["detail"]["decrypt"] = d["bit"], d
    row["status"] = ("OK_CERTIFIED" if r["certified"] else "UNCERTIFIED_BUDGET") if row["verified"] else "FAILED_VERIFICATION"
    return row


def print_lwe_result(res: Dict) -> None:
    print(f"\n== {res['name']}  {res['problem']}  k={res['k']} m={res['m']} q={res['q']}  -> {res['status']}")
    if res["answer"] is None:
        print("   ", res["note"])
        return
    d = res["detail"]
    print(f"   secreto s = {res['answer']}")
    print(f"   errores e = {d['errors_centered']}   (||e||^2 = {d['err_norm2']}, max|e| = {d['err_max']})")
    if res["bit"] is not None:
        c = d["decrypt"]
        print(f"   BIT CIFRADO = {res['bit']}   (d = v - <u,s> = {c['d']}; distancia a 0: {c['dist0']}, "
              f"a floor(q/2)={c['half']}: {c['dist1']})")
    print(f"   celosia incrustada dim={d['dim']}  backend={res['backend']}  t_LLL={res['t_reduce']:.4f}s  "
          f"t_enum={res['t_solve']:.4f}s  nodos={res['nodes']}  OPTIMO_CERTIFICADO={res['certified']}")
    print(f"   comprobacion: ecuaciones consistentes={d['consistent']}  CVP directo coincide={d.get('cvp_agrees')}  "
          f"hueco gauss/||e||={d['gap']:.2f}  log10(soluciones espurias esperadas)={d['spurious_log10']:.1f}")
    if res["note"]:
        print("   nota:", res["note"])
