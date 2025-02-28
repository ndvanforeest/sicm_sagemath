import numpy as np

load(
    "utils.sage",
    "utils1.4.sage",
    "utils1.5.sage",
    "utils1.6.sage",
)

def Lagrangian_to_energy(L):
    P = partial(L, 2)

    def f(local):
        v = velocity(local)
        return vector(P(local)) * v - L(local)

    return f


def Lagrangian_to_energy(L):
    P = partial(L, 2)
    return lambda local: Sum(Product(Compose(vector, P), velocity), Min(L))(
        local
    )

var("t x y", domain="real")

space = make_named_space(["r", "\\theta", "\\phi"])

def s_to_r(sperical_state):
    r, theta, phi = coordinate(sperical_state)[:]
    return vector(
        [r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)]
    )

def U(r):
    return function("V")(r)

def L_3D_central(m, V):
    def Lagrangian(local):
        return L_central_rectangular(m, V)(F_to_C(s_to_r)(local))

    return Lagrangian

show(partial(L_3D_central(m, U), 1)(space).simplify_full())

show(partial(L_3D_central(m, U), 2)(space).simplify_full())

def ang_mom_z(m):
    def f(rectangular_state):
        xyx = coordinate(rectangular_state)
        v = velocity(rectangular_state)
        return xyx.cross_product(m * v)[2]

    return f


show(Compose(ang_mom_z(m), F_to_C(s_to_r))(space).simplify_full())

show(Lagrangian_to_energy(L_central_spherical(m, U))(space).simplify_full())

var("G M0 M1 a")
assume(G > 0, M0 > 0, M1 > 0, a > 0)


def distance(x, y):
    return sqrt(square(x - y))


def angular_freq(M0, M1, a):
    return sqrt(G * (M0 + M1) / a ^ 3)


def V(a, M0, M1, m):
    Omega = angular_freq(M0, M1, a)
    a0, a1 = M1 / (M0 + M1) * a, M0 / (M0 + M1) * a

    def f(t, origin):
        pos0 = -a0 * vector([cos(Omega * t), sin(Omega * t)])
        pos1 = a1 * vector([cos(Omega * t), sin(Omega * t)])
        r0 = distance(origin, pos0)
        r1 = distance(origin, pos1)
        return -G * m * (M0 / r0 + M1 / r1)

    return f


def L0(m, V):
    def f(local):
        t, q, v = time(local), coordinate(local), velocity(local)
        return 1 / 2 * m * square(v) - V(t, q)

    return f

space = make_named_space(["x", "y"])
show(Lagrangian_to_energy(L0(m, V(a, M0, M1, m)))(space).simplify_full().expand())

def F_tilde(angle_x, angle_y, angle_z):
    def f(local):
        return (
            rotation_matrix([1, 0, 0], angle_x)
            * rotation_matrix([0, 1, 0], angle_y)
            * rotation_matrix([0, 0, 1], angle_z)
            * coordinate(local)
        )

    return f

space = make_named_space(["x", "y", "z"])

var("s t u")

def Rx(s):
    return F_tilde(s, 0, 0)(space)


show(Rx(s))
show(diff(Rx(s), s)(s=0))

def the_Noether_integral(local):
    L = L_central_rectangular(m, U)
    DF0 = jacobian(F_tilde(s, t, u)(local), (s, t, u))(s=0, t=0, u=0)
    return partial(L, 2)(local) * DF0


show(the_Noether_integral(space))
