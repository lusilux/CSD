"""
make_tables.py -- Genera las tablas de resultados de la memoria a partir de los CSV de solve.py.

Uso:   python make_tables.py tiempos.csv tiemposLwe.csv [--out tablas]
       (vale cualquier numero de CSV; las filas se clasifican solas: SVP/CVP -> celosias,
        LWE-* / RLWE-* -> LWE y Ring-LWE)

Salida: <out>.md   tablas en Markdown        <out>.tex   tablas en LaTeX (booktabs)
        <out>_anexo.md   respuestas completas (vectores, secretos, palabras de bits)
"""
import argparse
import csv
import re


# ------------------------------------------------------------------ utilidades
def natkey(name):
    m = re.match(r"Tier0*(\d+)_?T0*(\d+)", name)
    return (int(m.group(1)), int(m.group(2)), name) if m else (99, 0, name)


def f(row, key, default=None):
    v = row.get(key)
    return default if v in (None, "") else v


def num(row, key):
    v = f(row, key)
    return float(v) if v is not None else None


def fmt_t(t):
    if t is None:
        return "-"
    return f"{t * 1000:.2f} ms" if t < 1 else f"{t:.2f} s"


def tex_esc(s):
    for a, b in (("\\", r"\textbackslash{}"), ("_", r"\_"), ("%", r"\%"), ("&", r"\&"), ("#", r"\#")):
        s = s.replace(a, b)
    return s


class C:
    """Celda con version Markdown y version LaTeX."""
    def __init__(self, md, tex=None):
        self.md, self.tex = str(md), (tex if tex is not None else tex_esc(str(md)))


def render_md(headers, rows, align):
    sep = ["---:" if a == "r" else (":---:" if a == "c" else ":---") for a in align]
    out = ["| " + " | ".join(h.md for h in headers) + " |", "| " + " | ".join(sep) + " |"]
    out += ["| " + " | ".join(c.md for c in r) + " |" for r in rows]
    return "\n".join(out)


def render_tex(headers, rows, align, caption, label):
    out = [r"\begin{table}[ht]", r"\centering\small", f"\\caption{{{caption}}}", f"\\label{{{label}}}",
           r"\begin{tabular}{" + "".join(align) + "}", r"\toprule",
           " & ".join(h.tex for h in headers) + r" \\", r"\midrule"]
    out += [" & ".join(c.tex for c in r) + r" \\" for r in rows]
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(out)


def t_lll(r):
    return num(r, "t_lll") if f(r, "t_lll") else num(r, "t_reduce")      # compatible con CSV antiguos


def t_enum(r):
    return num(r, "t_enum") if f(r, "t_enum") else num(r, "t_solve")


# ------------------------------------------------------------------ tabla de celosias
def table_lattice(rows):
    H = [C("Ejercicio"), C("Tipo"), C("Dim."), C("Módulo"), C("Resultado", "Resultado"), C("Origen"),
         C("Swaps"), C("Nodos"), C("$t_{LLL}$", "$t_{LLL}$"), C("$t_{enum}$", "$t_{enum}$"), C("$t_{total}$", "$t_{total}$"),
         C("Cert.")]
    body, notes = [], []
    for r in sorted(rows, key=lambda r: natkey(r["name"])):
        svp = r["problem"] == "SVP"
        mod = f(r, "modulus", "N/A")
        degenerate = mod != "N/A"
        res = C(f"‖v‖² = {r['value']}" if svp else f"d² = {r['value']}",
                f"$\\|v\\|^2={r['value']}$" if svp else f"$d^2={r['value']}$")
        body.append([C(r["name"]), C(r["problem"]), C(r["dim"]),
                     C(f"{mod}†" if degenerate else "—", f"{mod}$^\\dagger$" if degenerate else "--"),
                     res, C(f(r, "source", "-")), C(f(r, "swaps", "-")), C(f(r, "nodes", "-")),
                     C(fmt_t(t_lll(r))), C(fmt_t(t_enum(r))), C(fmt_t(num(r, "t_total"))),
                     C("sí" if r["certified"] == "True" and r["status"] == "OK_CERTIFIED" else "NO")])
        if degenerate:
            notes.append(r["name"])
    foot = [("Origen = de dónde sale la respuesta (Babai, base reducida por LLL o enumeración). "
             "t_total incluye lectura del fichero y verificaciones independientes.",
             "Origen = de dónde sale la respuesta (Babai, base reducida por LLL o enumeración). "
             "$t_{total}$ incluye lectura del fichero y verificaciones independientes.")]
    if notes:
        foot.append(("† Instancia con módulo, resuelta con la lectura L(B)+nZ^k; la celosía resultante es degenerada "
                     "(determinante 1, es decir Z^k), por lo que el resultado es trivial: " + ", ".join(notes) + ".",
                     "$\\dagger$ Instancia con módulo, resuelta con la lectura $L(B)+n\\mathbb{Z}^k$; la celosía resultante es "
                     "degenerada (determinante 1, es decir $\\mathbb{Z}^k$), por lo que el resultado es trivial: "
                     + tex_esc(", ".join(notes)) + "."))
    return H, body, ["l", "l", "r", "r", "l", "l", "r", "r", "r", "r", "r", "c"], foot


# ------------------------------------------------------------------ tabla LWE / Ring-LWE
def table_lwe(rows):
    H = [C("Ejercicio"), C("Tarea"), C("$k$ ($n$)", "$k$ ($n$)"), C("$m$", "$m$"), C("$q$", "$q$"), C("Dim."),
         C("máx \\|e\\|", r"$\max|e|$"), C("Hueco"), C("Resultado"), C("Swaps"),
         C("$t_{LLL}$", "$t_{LLL}$"), C("$t_{enum}$", "$t_{enum}$"), C("$t_{total}$", "$t_{total}$"), C("Cert.")]
    names = {"LWE-rotura": "LWE rotura", "LWE-criptoanalisis": "LWE cripto.",
             "RLWE-rotura": "Ring-LWE rotura", "RLWE-criptoanalisis": "Ring-LWE cripto."}
    body = []
    for r in sorted(rows, key=lambda r: natkey(r["name"])):
        bit = f(r, "bit")
        if bit is None:
            res = "s (anexo)"
        elif len(bit) > 16:
            res = bit[:16] + "…"
        else:
            res = bit
        gap = num(r, "gap")
        body.append([C(r["name"]), C(names.get(r["problem"], r["problem"])), C(r["k"]), C(r["m"]), C(r["q"]),
                     C(r["dim"]), C(f(r, "e_max", "-")), C(f"{gap:.1f}" if gap is not None else "-"),
                     C(res, tex_esc(res).replace("…", r"\ldots")), C(f(r, "swaps", "-")),
                     C(fmt_t(t_lll(r))), C(fmt_t(t_enum(r))), C(fmt_t(num(r, "t_total"))),
                     C("sí" if r["status"] == "OK_CERTIFIED" else "NO")])
    foot = [("Hueco = λ1 esperado de la celosía / ‖e‖ (cuanto mayor, más seguro es que la solución sea única). "
             "Resultado: bit cifrado (LWE), palabra de bits con el bit i = coeficiente de x^i (Ring-LWE; completa en el "
             "anexo) o secreto s (rotura). En LWE, t_total incluye la verificación cruzada con CVP directo.",
             "Hueco = $\\lambda_1$ esperado de la celosía / $\\|e\\|$ (cuanto mayor, más seguro es que la solución sea única). "
             "Resultado: bit cifrado (LWE), palabra de bits con el bit $i$ = coeficiente de $x^i$ (Ring-LWE; completa en el "
             "anexo) o secreto $s$ (rotura). En LWE, $t_{total}$ incluye la verificación cruzada con CVP directo.")]
    return H, body, ["l", "l", "r", "r", "r", "r", "r", "r", "l", "r", "r", "r", "r", "c"], foot


# ------------------------------------------------------------------ anexo
def annex(rows):
    out = ["# Anexo: respuestas completas\n"]
    rows = sorted(rows, key=lambda r: (r["problem"] not in ("SVP", "CVP"), natkey(r["name"])))   # celosias primero
    for r in rows:
        line = f"- **{r['name']}** ({r['problem']}): "
        if f(r, "answer_text"):
            what = {"SVP": "vector más corto", "CVP": "punto más cercano"}.get(r["problem"], "secreto s")
            line += f"{what} = `{r['answer_text']}`"
        if f(r, "bit"):
            line += (f"; " if f(r, "answer_text") else "") + f"{'palabra de bits' if len(r['bit']) > 1 else 'bit cifrado'} = `{r['bit']}`"
        out.append(line)
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+")
    ap.add_argument("--out", default="tablas")
    a = ap.parse_args()
    rows = {}
    for path in a.csv:
        with open(path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rows[r["name"]] = r                    # si un nombre se repite, gana el ultimo CSV
    lat = [r for r in rows.values() if r["problem"] in ("SVP", "CVP")]
    lwe = [r for r in rows.values() if r["problem"].startswith(("LWE", "RLWE"))]
    other = len(rows) - len(lat) - len(lwe)
    md, tex = [], []
    for title, build, rs, cap, lab in (
            ("Celosías (SVP/CVP)", table_lattice, lat, "Resultados de los ejercicios de celosías.", "tab:celosias"),
            ("LWE y Ring-LWE", table_lwe, lwe, "Resultados de los ejercicios LWE y Ring-LWE.", "tab:lwe")):
        if not rs:
            continue
        H, body, align, foot = build(rs)
        md += [f"## {title}\n", render_md(H, body, align), "", *[f"*{x[0]}*  " for x in foot], ""]
        tex += [render_tex(H, body, align, cap + " " + " ".join(x[1] for x in foot), lab), ""]
    with open(a.out + ".md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))
    with open(a.out + ".tex", "w", encoding="utf-8") as fh:
        fh.write("% requiere \\usepackage{booktabs}\n" + "\n".join(tex))
    with open(a.out + "_anexo.md", "w", encoding="utf-8") as fh:
        fh.write(annex(lat + lwe))
    print(f"{len(lat)} celosías + {len(lwe)} LWE/Ring-LWE" + (f" ({other} filas ignoradas)" if other else "")
          + f" -> {a.out}.md, {a.out}.tex, {a.out}_anexo.md")


if __name__ == "__main__":
    main()
