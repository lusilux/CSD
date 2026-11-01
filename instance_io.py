"""Parser de los ficheros de instancia SVP/CVP (lineas '#' = comentario)."""
from typing import Dict, List


def parse_lattice_instance(path: str) -> Dict:
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = [l.strip() for l in f if l.strip() and not l.lstrip().startswith("#")]
    problem = lines[0].upper()                         # SVP | CVP
    if problem not in ("SVP", "CVP"):
        raise ValueError(f"{path}: la primera linea debe ser SVP o CVP (aqui: {lines[0]!r}); "
                         "los ficheros LWE todavia no estan soportados")
    target = [int(x) for x in lines[1].split()]        # "0 0" en SVP
    dim = int(lines[2])
    basis: List[List[int]] = [[int(x) for x in lines[3 + i].replace(",", " ").split()]
                              for i in range(dim)]
    modulus = lines[3 + dim]                           # numero o N/A
    solution = lines[4 + dim] if len(lines) > 4 + dim else "N/D"
    assert all(len(r) == dim for r in basis), "se esperaba base cuadrada dim x dim"
    return dict(problem=problem, target=target, dim=dim, basis=basis,
                modulus=None if modulus.upper() == "N/A" else int(modulus),
                solution=solution)
