import re
import tuples # see below

def show_expression(s):
    s = latex(s)
    s = re.sub(r"\\frac{\\partial}{\\partial t}", r"\\dot ", s)
    s = re.sub(r"\\left\(t\\right\)", r"", s)
    s = re.sub(
        r"\\frac\{\\partial\^\{2\}\}\{\(\\partial t\)\^\{2\}\}", r"\\ddot ", s
    )
    print(r"\[" + s + r"\]")


# def show_expression(s):
#     return r"\[" + latex(s) + r"\]"

def show_tuple(tup):
    res = r"\begin{align*}"
    for component in tup:
        res += "& " + latex(component) + r"\\"
    res += r"\end{align*}"
    return res

def show(s):
    if isinstance(s, tuples.Tuple):
        return show_tuple(s)
    return show_expression(s)
