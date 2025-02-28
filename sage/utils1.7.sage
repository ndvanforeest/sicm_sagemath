load(
    "utils.sage",
    "utils1.6.sage",
)

def Lagrangian_to_accelaration(L):
    def f(local):
        P = partial(L, 2)  # (local)
        F = partial(L, 1)  # (local)
        res = partial(P, 2)(local).solve_left(
            vector(F(local))
            - vector(partial(P, 0)(local))
            - partial(P, 1)(local) * velocity(local)
        )
        return res

    return f

def convert_to_expr(n):
    return SR(n)

def Lagrangian_to_state_derivative(L):
    accelaration = Lagrangian_to_acceleration(L)
    return lambda state: up(
        convert_to_expr(1), velocity(state), acceleration(state)
    )

def Lagrange_equations_first_order(L):
    def f(q, v):
        state_path = qv_to_state_path(q, v)(t)
        res = D(state_path)
        res -= Lagrangian_to_state_derivative(L)(state_path)
        return res

    return f

def state_to_list(state):
    return [time(state), *coordinate(state), *velocity(state)]

def evolve(state_derivative, ics, times):
    space = make_space("qq", dim=len(coordinate(ics)))
    soln = desolve_odeint(
        des=state_to_list(state_derivative(space)),
        ics=state_to_list(ics),
        times=times,
        dvars=state_to_list(space),
        rtol=1e-13,
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

var("m l g A omega")


def L_periodically_driven_pendulum(m, l, g, A, omega):
    ys = periodic_drive(A, omega, 0)

    def Lagrangian(local):
        return L_pend(m, l, g, ys)(local)

    return Lagrangian

def principal_value(cut_point):
    def f(x):
        return (x + cut_point) % (2 * np.pi) - cut_point

    return f
