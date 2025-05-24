load("utils1.5.sage")

var("t", domain="real")
var("g m", domain="positive")


def L_uniform_acceleration(m, g):
    def wrap_L_unif(local):
        x, y  = coordinate(local).list()
        v = velocity(local)
        T = 1 / 2 * m * square(v)
        V = m * g * y
        return T - V

    return wrap_L_unif

def L_central_rectangular(m, U):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        T = 1 / 2 * m * square(v)
        return T - U(sqrt(square(q)))

    return Lagrangian

def F_to_C(F):
    def wrap_F_to_C(local):
        return up(
            time(local),
            F(local),
            partial(F, 0)(local) + partial(F, 1)(local) * velocity(local),
        )

    return wrap_F_to_C

def p_to_r(local):
    r, phi = coordinate(local).list()
    return column_matrix([r * cos(phi), r * sin(phi)])

def L_central_polar(m, U):
    def Lagrangian(local):
        return compose(L_central_rectangular(m, U), F_to_C(p_to_r))(local)

    return Lagrangian

def L_free_rectangular(m):
    def Lagrangian(local):
        v = velocity(local)
        return 1 / 2 * m * square(v)

    return Lagrangian


def L_free_polar(m):
    def Lagrangian(local):
        return L_free_rectangular(m)(F_to_C(p_to_r)(local))

    return Lagrangian


def F(Omega):
    def f(local):
        t = time(local)
        r, theta = coordinate(local).list()
        return vector([r, theta + Omega * t])

    return f


def L_rotating_polar(m, Omega):
    def Lagrangian(local):
        return L_free_polar(m)(F_to_C(F(Omega))(local))

    return Lagrangian



def r_to_p(local):
    x, y = coordinate(local).list()
    return column_matrix([sqrt(x * x + y * y), atan(y / x)])


def L_rotating_rectangular(m, Omega):
    def Lagrangian(local):
        return L_rotating_polar(m, Omega)(F_to_C(r_to_p)(local))

    return Lagrangian

def dp_coordinates(l, ys):
    "From theta to x, y coordinates."
    def wrap_dp(local):
        t = time(local)
        theta = coordinate(local)[0, 0]
        return column_matrix([l * sin(theta), ys(t) - l * cos(theta)])

    return wrap_dp

def L_pend(m, l, g, ys):
    def wrap_L_pend(local):
        return L_uniform_acceleration(m, g)(
            F_to_C(dp_coordinates(l, ys))(local)
        )

    return wrap_L_pend
