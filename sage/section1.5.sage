import numpy as np

load(
    "utils.sage",
    "utils1.4.sage",
    "utils1.5.sage",
)

var("t x y", domain="real")

var('k m')
assume(k > 0, m > 0)
space = make_named_space(["q_x"])

show(L_harmonic(m, k)(space))

show(partial(L_harmonic(m, k), 1)(space))

show(partial(L_harmonic(m, k), 2)(space))

q = vector([literal_function("q_x"), literal_function("q_y")])

show(partial(L_harmonic(m, k), 1)(Gamma(q)(t)))

show(Compose(partial(L_harmonic(m, k), 1), Gamma(q))(t))

show(D(partial(L_harmonic(m, k), 2)(Gamma(q)(t))))

space = make_named_space(["\\xi", "\\eta"])
q = vector([literal_function("\\xi"), literal_function("\\eta")])

var("mu", domain="positive")

def L_orbital(m, mu):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        return (1 / 2) * m * square(v) + mu / sqrt(square(q))

    return Lagrangian

L = L_orbital(m, mu)
show(L(space))

show(partial(L, 1)(Gamma(q)(t)))

show(partial(L, 2)(Gamma(q)(t)))

space = make_named_space(["\\theta"])
q = vector([literal_function("\\theta")])

L = L_planar_pendulum(m, g, l)
show(L(space))

show(partial(L, 1)(Gamma(q)(t)))

show(partial(L, 2)(Gamma(q)(t)))

L = L_Henon_Heiles(m)
space = make_space("x", dim=2)
show(L(space))

show(partial(L, 1)(space))

show(partial(L, 2)(space))

var('R', domain="positive")


def L_sphere(m, R):
    def Lagrangian(local):
        q = coordinate(local)
        theta, phi = q[:]
        v = velocity(local)
        alpha, beta = v[:]
        L = m * R * (square(alpha) + square(beta * sin(theta))) / 2
        return L

    return Lagrangian

space = make_named_space(["\\phi", "\\theta"])
L = L_sphere(m, R)
show(L(space))

show(partial(L, 1)(space))

show(partial(L, 2)(space))

space = make_space("x", dim=3)
var("a b c a0 b0 c0", domain="real")
test_path = vector([a * t + a0, b * t + b0, c * t + c0])

l_eq = Lagrange_equations(L_free_particle(m))(test_path)
show(l_eq)
show(l_eq(t))

l_eq = Lagrange_equations(L_free_particle(m))(
    [literal_function("x"), literal_function("y"), literal_function('z')]
)
show(l_eq(t))

var("A phi omega", domain="real")
assume(A > 0)

proposed_path = vector([A * cos(omega * t + phi)])

l_eq = Lagrange_equations(L_harmonic(m, k))(proposed_path)(t)
show(l_eq)

show(l_eq[0][0])

show(l_eq[0][0].factor())

space = make_named_space(["r", "\\phi"])
var("G m m1 m2")
assume(G > 0, m > 0, m1 > 0, m2 > 0)


def L_central_polar(m, V):
    def Lagrangian(local):
        r, phi = coordinate(local)
        # r, phi = q[:]
        qdot = velocity(local)
        rdot, phidot = velocity(local)  # qdot[:]
        T = 1 / 2 * m * (square(rdot) + square(r * phidot))
        return T - V(r)

    return Lagrangian


def gravitational_energy(G, m1, m2):
    def f(r):
        return -G * m1 * m2 / r

    return f

V = gravitational_energy(G, m1, m2)
L = L_central_polar(m, gravitational_energy(G, m1, m2))
show(L(space))

l_eq = Lagrange_equations(L_central_polar(m, gravitational_energy(G, m1, m2)))(
    [literal_function("r"), literal_function("\\phi")]
)(t)

show(l_eq[0][1].factor() == 0)

show(l_eq[0][0] == 0)
