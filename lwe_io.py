"""
lwe_io.py -- Lector de los ficheros LWE de la asignatura (lineas '#' = comentario).

Formato (tras quitar comentarios y lineas vacias):
    1  tarea:  'criptoanalisis' (descifrar un bit) | 'rotura' (recuperar el vector de secretos)
    2  criptograma:  {u1,...,uk}, v      (N/A en 'rotura')
    3  k  = numero de incognitas (el enunciado lo llama 'dimension del sistema')
    4  q  = modulo
    5  m lineas {a_i1,...,a_ik}          (coeficientes de las m ecuaciones)
    6  1 linea  {b_1,...,b_m}            (terminos independientes)
    7  solucion (o N/D)
Modelo (PDF, Regev): b_i = <a_i, s> + e_i (mod q), con errores e_i pequenos.
"""
import re
from typing import Dict, List

TASKS = {"criptoanalisis": "criptoanalisis", "criptoanálisis": "criptoanalisis", "rotura": "rotura"}


def _data_lines(path: str) -> List[str]:
    with open(path, encoding="utf-8", errors="replace") as f:
        return [l.strip() for l in f if l.strip() and not l.lstrip().startswith("#")]


def is_lwe_file(path: str) -> bool:
    """True si la primera linea de datos es una tarea LWE (y no SVP/CVP)."""
    lines = _data_lines(path)
    return bool(lines) and lines[0].strip().lower() in TASKS


def _ints(s: str) -> List[int]:
    return [int(x) for x in re.findall(r"-?\d+", s)]


def parse_lwe_instance(path: str) -> Dict:
    lines = _data_lines(path)
    task = TASKS.get(lines[0].strip().lower())
    if task is None:
        raise ValueError(f"{path}: tarea LWE desconocida: {lines[0]!r}")
    cipher = None
    if lines[1].strip().upper() != "N/A":
        nums = _ints(lines[1])
        cipher = (nums[:-1], nums[-1])
    k, q = int(lines[2]), int(lines[3])
    brace_idx = [i for i in range(4, len(lines)) if lines[i].startswith("{")]
    if len(brace_idx) < 2:
        raise ValueError(f"{path}: faltan las ecuaciones o los terminos independientes")
    rows = [_ints(lines[i]) for i in brace_idx[:-1]]
    b = _ints(lines[brace_idx[-1]])
    solution = lines[brace_idx[-1] + 1] if brace_idx[-1] + 1 < len(lines) else "N/D"
    if any(len(r) != k for r in rows):
        raise ValueError(f"{path}: todas las ecuaciones deben tener {k} coeficientes")
    if len(b) != len(rows):
        raise ValueError(f"{path}: {len(rows)} ecuaciones pero {len(b)} terminos independientes")
    if cipher is not None and len(cipher[0]) != k:
        raise ValueError(f"{path}: el criptograma debe tener {k} coeficientes")
    return dict(problem="LWE", task=task, cipher=cipher, k=k, q=q, A=rows, b=b, m=len(rows),
                solution=solution)
