import numpy as np

load(
    "utils.sage",
    "utils1.4.sage",
    "utils1.5.sage",
    "utils1.6.sage",
)

var("t x y", domain="real")

space = make_named_space(["x", "y"])

l_eq = Lagrange_equations(L_uniform_acceleration(m, g))(
    [literal_function("x"), literal_function("y")]
)
show(l_eq)

def U(r):
    return 1 / r

show(
    Lagrange_equations(L_central_rectangular(m, U))(
        [literal_function("x"), literal_function("y")]
    )
)

def U(r):
    return function("V")(r)

show(
    Lagrange_equations(L_central_rectangular(m, U))(
        [literal_function("x"), literal_function("y")]
    )
)

space = make_named_space(["r", "\\phi"])
show(p_to_r(space))

show((partial(p_to_r, 0)(space)))

show((partial(p_to_r, 1)(space)))

show(F_to_C(p_to_r)(space))

show(L_central_polar(m, U)(space).simplify_full())

expr = Lagrange_equations(L_central_polar(m, U))(
    [literal_function("r"), literal_function("\\phi")]
).simplify_full().expand()

show(expr[0][0])
show(expr[0][1])

space = make_named_space(["x", "y"])
var("m Omega r")
expr = L_rotating_rectangular(m, Omega)(space).simplify_full()

show(expr)

space = make_named_space(["\\theta"])
ys = literal_function("y")

expr = L_pend(m, l, g, ys)(space).simplify_full()
show(expr)
