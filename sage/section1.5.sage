load("utils1.5.sage")

t = var("t", domain="real")

load("utils1.4.sage")
k, m = var('k m', domain="positive")
q = column_path([literal_function("x")])

L = L_harmonic(m, k)
show(L(Gamma(q)(t)))

show(partial(L, 1)(Gamma(q)(t)))

show(partial(L, 2)(Gamma(q)(t)))

show(compose(partial(L, 1), Gamma(q))(t))
show(compose(partial(L, 2), Gamma(q))(t))

show(D(compose(partial(L, 2), Gamma(q)))(t))

q = column_path([literal_function("xi"), literal_function("eta")])

var("mu", domain="positive")

def L_orbital(m, mu):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        return (1 / 2) * m * square(v) + mu / sqrt(square(q))

    return Lagrangian

L = L_orbital(m, mu)
show(L(Gamma(q)(t)))

show(partial(L, 1)(Gamma(q)(t)))

show(partial(L, 2)(Gamma(q)(t)))

q = column_path([literal_function("theta")])

L = L_planar_pendulum(m, g, l)
show(L(Gamma(q)(t)))

show(partial(L, 1)(Gamma(q)(t)))

show(partial(L, 2)(Gamma(q)(t)))

L = L_Henon_Heiles(m)
q = column_path([literal_function("x"), literal_function("y")])
show(L(Gamma(q)(t)))

show(partial(L, 1)(Gamma(q)(t)))

show(partial(L, 2)(Gamma(q)(t)))

var('R', domain="positive")


def L_sphere(m, R):
    def Lagrangian(local):
        theta, phi = coordinate(local).list()
        alpha, beta = velocity(local).list()
        L = m * R * (square(alpha) + square(beta * sin(theta))) / 2
        return L

    return Lagrangian

q = column_path([literal_function("phi"), literal_function("theta")])
L = L_sphere(m, R)

show(L(Gamma(q)(t)))

show(partial(L, 1)(Gamma(q)(t)))

show(partial(L, 2)(Gamma(q)(t)))

q = column_path(
    [
        literal_function("x"),
        literal_function("y"),
    ]
)

L = L_free_particle(m)
show(compose(partial(L, 1), Gamma(q))(t))
show(compose(partial(L, 2), Gamma(q))(t))
show(D(compose(partial(L, 2), Gamma(q)))(t))
show(
    (D(compose(partial(L, 2), Gamma(q))) - compose(partial(L, 1), Gamma(q)))(t)
)

var("a b c a0 b0 c0", domain="real")
test_path = lambda t: column_matrix([a * t + a0, b * t + b0, c * t + c0])

l_eq = Lagrange_equations(L_free_particle(m))(test_path)
show(l_eq(t))

q = column_path([literal_function("x")])
l_eq = Lagrange_equations(L_free_particle(m))(q)
show(l_eq(t))

var("A phi omega", domain="real")
assume(A > 0)

proposed_path = lambda t: vector([A * cos(omega * t + phi)])

l_eq = Lagrange_equations(L_harmonic(m, k))(proposed_path)(t)
show(l_eq)

show(l_eq[0][0])

show(l_eq[0, 0].factor())

var("G m m1 m2", domain="positive")


def L_central_polar(m, V):
    def Lagrangian(local):
        r, phi = coordinate(local).list()
        rdot, phidot = velocity(local).list()
        T = 1 / 2 * m * (square(rdot) + square(r * phidot))
        return T - V(r)

    return Lagrangian


def gravitational_energy(G, m1, m2):
    def f(r):
        return -G * m1 * m2 / r

    return f

q = column_path([literal_function("r"), literal_function("phi")])
V = gravitational_energy(G, m1, m2)
L = L_central_polar(m, V)
show(L(Gamma(q)(t)))

l_eq = Lagrange_equations(L)(q)(t)

show(l_eq[0, 1] == 0)

show(l_eq[0, 0] == 0)
