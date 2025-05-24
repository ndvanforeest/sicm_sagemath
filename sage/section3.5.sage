import numpy as np

load("utils1.7.sage", "utils3.1.sage")

var("t", domain="real")

q = column_path([literal_function("theta")])
p = row_path([literal_function(r"p_theta")])
H_state = qp_to_H_state_path(q, p)(t)
show(H_state)

_ = var("A g l m omega", domain="positive")


H = Lagrangian_to_Hamiltonian(
    L_periodically_driven_pendulum(m, l, g, A, omega)
)


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
p.save(f'piet.png')

quit()
p.save(f'./../figures/hamiltonian_driven_pendulum_0.01.png')

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
p.save(f'./../figures/hamiltonian_driven_pendulum_0.005.png')
