var("m g l")
assume(m > 0, g > 0, l > 0)


def L_planar_pendulum(m, g, l):
    def Lagrangian(local):
        theta = coordinate(local)[0]
        theta_dot = velocity(local)
        T = (1 / 2) * m * l ^ 2 * square(theta_dot)
        V = m * g * l * (1 - cos(theta))
        return T - V

    return Lagrangian

def L_Henon_Heiles(m):
    def Lagrangian(local):
        q = coordinate(local)
        x, y = q[:]
        v = velocity(local)
        T = (1 / 2) * square(v)
        V = 1 / 2 * (square(x) + square(y)) + square(x) * y - y**3 / 3
        return T - V

    return Lagrangian

def Lagrange_equations(Lagrangian):
    def f(q):
        return Sum(
            Compose(D, partial(Lagrangian, 2), Gamma(q)),
            Min(Compose(partial(Lagrangian, 1), Gamma(q))),
        )

    return f
