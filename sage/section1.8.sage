load("utils1.8.sage")

var("t", domain=RR)

def Lagrangian_to_energy(L):
    P = partial(L, 2)
    LL = Function(lambda local: L(local))
    return lambda local: (P * velocity - LL)(local)

q = column_path(
    [
        literal_function("r"),
        literal_function("theta"),
        literal_function("phi"),
    ]
)

def s_to_r(sperical_state):
    r, theta, phi = coordinate(spherical_state).list()
    return vector(
        [r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)]
    )

show(velocity(F_to_C(s_to_r)(Gamma(q)(t))).simplify_full())

V = Function(lambda r: function("V")(r))

def L_3D_central(m, V):
    def Lagrangian(local):
        return L_central_rectangular(m, V)(F_to_C(s_to_r)(local))

    return Lagrangian

show(partial(L_3D_central(m, V), 1)(Gamma(q)(t)).simplify_full())

show(partial(L_3D_central(m, V), 2)(Gamma(q)(t)).simplify_full())

def ang_mom_z(m):
    def f(rectangular_state):
        xyx = vector(coordinate(rectangular_state))
        v = vector(velocity(rectangular_state))
        return xyx.cross_product(m * v)[2]

    return f


show(compose(ang_mom_z(m), F_to_C(s_to_r))(Gamma(q)(t)).simplify_full())

show(Lagrangian_to_energy(L_3D_central(m, V))(Gamma(q)(t)).simplify_full())

var("G M0 M1 a", domain="positive")


def distance(x, y):
    return sqrt(square(x - y))


def angular_freq(M0, M1, a):
    return sqrt(G * (M0 + M1) / a ^ 3)


def V(a, M0, M1, m):
    Omega = angular_freq(M0, M1, a)
    a0, a1 = M1 / (M0 + M1) * a, M0 / (M0 + M1) * a

    def f(t, origin):
        pos0 = -a0 * column_matrix([cos(Omega * t), sin(Omega * t)])
        pos1 = a1 * column_matrix([cos(Omega * t), sin(Omega * t)])
        r0 = distance(origin, pos0)
        r1 = distance(origin, pos1)
        return -G * m * (M0 / r0 + M1 / r1)

    return f

def L0(m, V):
    def f(local):
        t, q, v = time(local), coordinate(local), velocity(local)
        return 1 / 2 * m * square(v) - V(t, q)

    return f

q = column_path([literal_function("x"), literal_function("y")])
expr = (sqrt(G*M0 + G*M1)*t) / a^(3/2)
A = var('A')

show(
    Lagrangian_to_energy(L0(m, V(a, M0, M1, m)))(Gamma(q)(t))
    .simplify_full()
    .expand()
    .subs({expr: A})
)

def F_tilde(angle_x, angle_y, angle_z):
    def f(local):
        return (
            rotation_matrix([1, 0, 0], angle_x)
            * rotation_matrix([0, 1, 0], angle_y)
            * rotation_matrix([0, 0, 1], angle_z)
            * coordinate(local)
        )

    return f

q = column_path(
    [literal_function("x"), literal_function("y"), literal_function("z")]
)

def Rx(s):
    return lambda local: F_tilde(s, 0, 0)(local)


s, u, v = var("s u v")
latex.matrix_delimiters(left='[', right=']')
latex.matrix_column_alignment("c")
show(Rx(s)(Gamma(q)(t)))
show(diff(Rx(s)(Gamma(q)(t)), s)(s=0))

U = Function(lambda r: function("U")(r))


def the_Noether_integral(local):
    L = L_central_rectangular(m, U)
    Ftilde = lambda x: F_tilde(*x)(local)
    DF0 = Jacobian(Ftilde)([s, u, v], [s, u, v])(s=0, u=0, v=0)
    return partial(L, 2)(local) * DF0

show(the_Noether_integral(Gamma(q)(t)).simplify_full())
