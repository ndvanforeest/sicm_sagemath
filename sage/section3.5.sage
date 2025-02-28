import numpy as np

load(
    "utils.sage",
    "utils1.7.sage",
    "utils3.1.sage",
)

var("t x y", domain="real")

q = vector([literal_function("\\theta")])
p = vector([literal_function(r"p_\theta")])
H_state = qp_to_H_state_path(q, p)(t)
show(H_state)

H = Lagrangian_to_Hamiltonian(
    L_periodically_driven_pendulum(m, l, g, A, omega)
)
expr = H(H_state).simplify_full()
show(expr.expand())

DH = Hamiltonian_to_state_derivative(H)(H_state)
show(DH[1].simplify_full())
show(DH[2].simplify_full())

def H_pend_sysder(m, l, g, A, omega):
    Hamiltonian = Lagrangian_to_Hamiltonian(
        L_periodically_driven_pendulum(m, l, g, A, omega)
    )

    def f(state):
        return Hamiltonian_to_state_derivative(Hamiltonian)(state)

    return f

times = srange(0, 100, 0.01, include_endpoint=True)
soln = evolve(
    H_pend_sysder(m=1, l=1, g=9.8, A=0.1, omega=2 * sqrt(9.8)),
    ics=up(0, vector([1]), vector([0])),
    times=times,
)
thetas = principal_value(np.pi)(soln[:, 1])
thetadots = soln[:, 2]
pp = list(zip(thetas, thetadots))
p = points(pp, color='blue', size=3)
p.save(f'./figures/hamiltonian_driven_pendulum_0.01.png')

times = srange(0, 100, 0.005, include_endpoint=True)
soln = evolve(
    H_pend_sysder(m=1, l=1, g=9.8, A=0.1, omega=2 * sqrt(9.8)),
    ics=up(0, vector([1]), vector([0])),
    times=times,
)
thetas = principal_value(np.pi)(soln[:, 1])
thetadots = soln[:, 2]
pp = list(zip(thetas, thetadots))
p = points(pp, color='blue', size=3)
p.save(f'./figures/hamiltonian_driven_pendulum_0.005.png')
