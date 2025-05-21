load("utils1.6.sage")

def Lagrangian_to_acceleration(L):
    def f(local):
        P = partial(L, 2)
        F = compose(transpose, partial(L, 1))
        M = (F - partial(P, 0)) - partial(P, 1) * velocity
        return partial(P, 2)(local).solve_right(M(local))

    return f

def convert_to_expr(n):
    return SR(n)

def Lagrangian_to_state_derivative(L):
    acceleration = Lagrangian_to_acceleration(L)
    return lambda state: up(
        convert_to_expr(1), velocity(state), acceleration(state)
    )

def qv_to_state_path(q, v):
    return lambda t: up(t, q(t), v(t))

def Lagrange_equations_first_order(L):
    def f(q, v):
        state_path = qv_to_state_path(q, v)
        res = D(state_path)
        res -= compose(Lagrangian_to_state_derivative(L), state_path)
        return res

    return f

def make_dummy_vector(name, dim):
    return column_matrix([var(f"{name}{i}", domain=RR) for i in range(dim)])

def evolve(state_derivative, ics, times):
    dim = coordinate(ics).nrows()
    coordinates = make_dummy_vector("q", dim)
    velocities = make_dummy_vector("v", dim)
    space = up(t, coordinates, velocities)
    soln = desolve_odeint(
        des=state_derivative(space).list(),
        ics=ics.list(),
        times=times,
        dvars=space.list(),
        atol=1e-13,
    )
    return soln

def state_advancer(state_derivative, ics, T):
    init_time = time(ics)
    times = [init_time, init_time + T]
    soln = evolve(state_derivative, ics, times)
    return soln[-1]

def periodic_drive(amplitude, frequency, phase):
    def f(t):
        return amplitude * cos(frequency * t + phase)

    return f

_ = var("m l g A omega")


def L_periodically_driven_pendulum(m, l, g, A, omega):
    ys = periodic_drive(A, omega, 0)

    def Lagrangian(local):
        return L_pend(m, l, g, ys)(local)

    return Lagrangian

def principal_value(cut_point):
    def f(x):
        return (x + cut_point) % (2 * np.pi) - cut_point

    return f
