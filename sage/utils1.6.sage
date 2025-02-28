var("g m")

def L_uniform_acceleration(m, g):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        y = q[1]
        T = 1 / 2 * m * v * v
        V = m * g * y
        return T - V

    return Lagrangian

def L_central_rectangular(m, U):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        T = 1 / 2 * m * square(v)
        return T - U(sqrt(square(q)))

    return Lagrangian

def F_to_C(F):
    def f(local):
        return up(
            time(local),
            F(local),
            vector(partial(F, 0)(local))
            + partial(F, 1)(local) * velocity(local),
        )

    return f

def p_to_r(local):
    r, phi = coordinate(local)
    return vector([r * cos(phi), r * sin(phi)])

def L_central_polar(m, U):
    def Lagrangian(local):
        return Compose(L_central_rectangular(m, U), F_to_C(p_to_r))(local)
        # return L_central_rectangular(m, U)(F_to_C(p_to_r)(local))

    return Lagrangian

def L_free_rectangular(m):
    def Lagrangian(local):
        v = velocity(local)
        return 1 / 2 * m * v * v

    return Lagrangian


def L_free_polar(m):
    def Lagrangian(local):
        return L_free_rectangular(m)(F_to_C(p_to_r)(local))

    return Lagrangian


def F(Omega):
    def f(local):
        t = time(local)
        r, theta = coordinate(local)
        return vector([r, theta + Omega * t])

    return f


def L_rotating_polar(m, Omega):
    def Lagrangian(local):
        return L_free_polar(m)(F_to_C(F(Omega))(local))

    return Lagrangian



# atan2(y/x) is not accepted when computing the L-ea
def r_to_p(local):
    x, y = coordinate(local)
    return vector([sqrt(x * x + y * y), atan(y / x)])


def L_rotating_rectangular(m, Omega):
    def Lagrangian(local):
        return L_rotating_polar(m, Omega)(F_to_C(r_to_p)(local))

    return Lagrangian

def dp_coordinates(l, ys):
    "From theta to x, y coordinates."
    def f(local):
        t = time(local)
        theta = coordinate(local)[0]
        return vector([l * sin(theta), ys(t=t) - l * cos(theta)])

    return f

var('g l m')
def L_pend(m, l, g, ys):
    def Lagrangian(local):
        return L_uniform_acceleration(m, g)(
            F_to_C(dp_coordinates(l, ys))(local)
        )

    return Lagrangian
