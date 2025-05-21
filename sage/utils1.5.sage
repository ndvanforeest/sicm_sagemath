load("utils1.4.sage")

var("m g l", domain="positive")


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
        x, y = coordinate(local).list()
        v = velocity(local)
        T = (1 / 2) * square(v)
        V = 1 / 2 * (square(x) + square(y)) + square(x) * y - y**3 / 3
        return T - V

    return Lagrangian

def Lagrange_equations(L):
    def f(q):
        return D(compose(partial(L, 2), Gamma(q))) - compose(
            partial(L, 1), Gamma(q)
        )

    return f
