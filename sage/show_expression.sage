import re

latex.matrix_delimiters(left='[', right=']')
latex.matrix_column_alignment("c")

def simplify_latex(s):
    s = re.sub(r"\\frac{\\partial}{\\partial t}", r"\\dot ", s)
    s = re.sub(r"\\left\(t\\right\)", r"", s)
    s = re.sub(
        r"\\frac\{\\partial\^\{2\}\}\{\(\\partial t\)\^\{2\}\}",
        r"\\ddot ",
        s,
    )
    return s

def show_expression(s, simplify=True):
    s = latex(s)
    if simplify:
        s = simplify_latex(s)
    res = r"\begin{dmath*}"
    res += "\n" + s + "\n"
    res += r"\end{dmath*}"
    print(res)

def show(s, simplify=True):
    return show_expression(s, simplify)
