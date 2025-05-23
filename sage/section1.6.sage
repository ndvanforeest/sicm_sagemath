load("utils1.6.sage")

q = column_path([literal_function("x"), literal_function("y")])
l_eq = Lagrange_equations(L_uniform_acceleration(m, g))(q)
show(l_eq(t))

def U(r):
    return 1 / r

show(Lagrange_equations(L_central_rectangular(m, U))(q)(t))

U = Function(lambda x: function("U")(x))
show(Lagrange_equations(L_central_rectangular(m, U))(q)(t))

r = literal_function("r")
phi = literal_function("phi")
q = column_path([r, phi])
show(p_to_r(Gamma(q)(t)))

show((partial(p_to_r, 0)(Gamma(q)(t))))

show((partial(p_to_r, 1)(Gamma(q)(t))))

show(F_to_C(p_to_r)(Gamma(q)(t)))

# show(L_central_polar(m, U)(Gamma(q)(t)))
show(L_central_polar(m, U)(Gamma(q)(t)).simplify_full())

expr = Lagrange_equations(L_central_polar(m, U))(q)(t)
show(expr.simplify_full().expand())

_ = var("Omega", domain="positive")
q_xy = column_path([literal_function("x"), literal_function("y")])
expr = L_rotating_rectangular(m, Omega)(Gamma(q_xy)(t)).simplify_full()

show(expr)

_ = var("l", domain="positive")

theta = column_path([literal_function("theta")])
ys = literal_function("y")

expr = L_pend(m, l, g, ys)(Gamma(theta)(t)).simplify_full()
show(expr)
