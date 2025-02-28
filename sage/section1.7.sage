import numpy as np

load(
    "utils.sage",
    "utils1.4.sage",
    "utils1.5.sage",
    "utils1.7.sage",
)

var("t x y", domain="real")

var("m k")
space = make_space("x", 2)

L = L_harmonic(m, k)

show(Lagrangian_to_accelaration(L)(space))

show_tuple(Lagrangian_to_state_derivative(L)(space))

def harmonic_state_derivative(m, k):
    return Lagrangian_to_state_derivative(L_harmonic(m, k))

show(harmonic_state_derivative(m, k)(space))

res = Lagrange_equations_first_order(L)(
    vector([literal_function("x"), literal_function("y")]),
    vector([literal_function("v_x"), literal_function("v_y")]),
)
show_tuple(res)

state_advancer(
    harmonic_state_derivative(m=2, k=1),
    ics=up(1, vector([1, 2]), vector([3, 4])),
    T=10,
)

space = make_named_space(["\\theta"])
show(L_periodically_driven_pendulum(m, l, g, A, omega)(space).simplify_full())

expr = Lagrange_equations(L_periodically_driven_pendulum(m, l, g, A, omega))(
    [literal_function("\\theta")]
).simplify_full()
show(expr)

show(
    Lagrangian_to_accelaration(
        L_periodically_driven_pendulum(m, l, g, A, omega)
    )(space).simplify_full()
)

def pend_state_derivative(m, l, g, A, omega):
    return Lagrangian_to_state_derivative(
        L_periodically_driven_pendulum(m, l, g, A, omega)
    )

expr = pend_state_derivative(m, l, g, A, omega)(space)
show(velocity(expr).simplify_full())

def plot_driven_pendulum(A, T, step_size=0.01):
    times = srange(0, T, step_size, include_endpoint=True)
    soln = evolve(
        pend_state_derivative(m=1, l=1, g=9.8, A=A, omega=2 * sqrt(9.8)),
        ics=up(0, vector([1]), vector([0])),
        times=times,
    )
    thetas = soln[:, 1]
    pp = list(zip(times, thetas))
    p = points(pp, color='blue', size=3)
    p.save(f'figures/driven_pendulum_{A:.2f}.png')

    thetas = principal_value(np.pi)(thetas)
    pp = list(zip(times, thetas))
    p = points(pp, color='blue', size=3)
    p.save(f'figures/driven_pendulum_{A:.2f}_principal_value.png')

    thetadots = soln[:, 2]
    pp = list(zip(thetas, thetadots))
    p = points(pp, color='blue', size=3)
    p.save(f'figures/driven_pendulum_{A:.2f}_trajectory.png')

plot_driven_pendulum(A=0.1, T=100, step_size=0.01)
