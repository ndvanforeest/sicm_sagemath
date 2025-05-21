load("utils1.7.sage")

var("t", domain=RR)

q = path_function([literal_function("x"), literal_function("y")])
local = Gamma(q)(t)
m, k = var("m k", domain="positive")
L = L_harmonic(m, k)
show(L(local))

F = compose(transpose, partial(L, 1))
show(F(local))
P = partial(L, 2)
show((F - partial(P, 0))(local))

show((partial(P, 1) * velocity)(local))

show((F - partial(P, 0) - partial(P, 1) * velocity)(local))

show(Lagrangian_to_acceleration(L)(local))

show(Lagrangian_to_state_derivative(L)(local))

def harmonic_state_derivative(m, k):
    return Lagrangian_to_state_derivative(L_harmonic(m, k))

show(harmonic_state_derivative(m, k)(local))

res = Lagrange_equations_first_order(L_harmonic(m, k))(
    path_function([literal_function("x"), literal_function("y")]),
    path_function([literal_function("v_x"), literal_function("v_y")]),
)
show(res(t))

state_advancer(
    harmonic_state_derivative(m=2, k=1),
    ics=up(0, column_matrix([1, 2]), column_matrix([3, 4])),
    T=10,
)

q = path_function([literal_function("theta")])
show(
    L_periodically_driven_pendulum(m, l, g, A, omega)(
        Gamma(q)(t)
    ).simplify_full()
)

expr = Lagrange_equations(L_periodically_driven_pendulum(m, l, g, A, omega))(
    q
)(t).simplify_full()
show(expr)

show(
    Lagrangian_to_acceleration(
        L_periodically_driven_pendulum(m, l, g, A, omega)
    )(Gamma(q)(t)).simplify_full()
)

def pend_state_derivative(m, l, g, A, omega):
    return Lagrangian_to_state_derivative(
        L_periodically_driven_pendulum(m, l, g, A, omega)
    )

expr = pend_state_derivative(m, l, g, A, omega)(Gamma(q)(t))
show(time(expr))
show(coordinate(expr).simplify_full())
show(velocity(expr).simplify_full())

def plot_driven_pendulum(A, T, step_size=0.01):
    times = srange(0, T, step_size, include_endpoint=True)
    soln = evolve(
        pend_state_derivative(m=1, l=1, g=9.8, A=A, omega=2 * sqrt(9.8)),
        ics=up(0, column_matrix([1]), column_matrix([0])),
        times=times,
    )
    thetas = soln[:, 1]
    pp = list(zip(times, thetas))
    p = points(pp, color='blue', size=3)
    p.save(f'../figures/driven_pendulum_{A:.2f}.png')

    thetas = principal_value(np.pi)(thetas)
    pp = list(zip(times, thetas))
    p = points(pp, color='blue', size=3)
    p.save(f'../figures/driven_pendulum_{A:.2f}_principal_value.png')

    thetadots = soln[:, 2]
    pp = list(zip(thetas, thetadots))
    p = points(pp, color='blue', size=3)
    p.save(f'../figures/driven_pendulum_{A:.2f}_trajectory.png')

plot_driven_pendulum(A=0.1, T=100, step_size=0.005)
