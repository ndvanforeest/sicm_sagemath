import numpy as np


load(
    "utils.sage",
    "utils1.5.sage",
    "utils1.6.sage",
    "utils3.1.sage",
)

var("t x y", domain="real")

def V(q):
    return function("U")(*q)

var("m")
q = vector([literal_function("q_x"), literal_function("q_y")])
p = vector([literal_function("p_x"), literal_function("p_y")])

H_state = qp_to_H_state_path(q, p)(t)
show(H_state)

H = H_rectangular
show(H(m, V)(H_state))

show(partial(H(m, V), 1)(H_state))

show_tuple(Hamiltonian_to_state_derivative(H(m, V))(H_state))

show_tuple(Hamilton_equations(H(m, V))(q, p))

res = Lagrangian_to_Hamiltonian(L_central_rectangular(m, V))(H_state)
show(res)

show(res.simplify_full())

var("m g l")
q = vector([literal_function("q_x")])
p = vector([literal_function("p_x")])

# space = make_named_space(["\\theta"])
H_state = qp_to_H_state_path(q, p)(t)
show(Lagrangian_to_Hamiltonian(L_planar_pendulum(m, g, l))(H_state))

q = vector([literal_function("q_x"), literal_function("q_y")])
p = vector([literal_function("p_x"), literal_function("p_y")])
H_state = qp_to_H_state_path(q, p)(t)
show(Lagrangian_to_Hamiltonian(L_Henon_Heiles(m))(H_state))

def L_sphere(m, R):
    def Lagrangian(local):
        theta, phi = coordinate(local)
        thetadot, phidot = velocity(local)
        return 1 / 2 * m * R ^ 2 * (
            square(thetadot) + square(phidot * sin(theta))
        )

    return Lagrangian


var("R")

space = make_named_space(["\\theta", "\\phi"])
show(Lagrangian_to_Hamiltonian(L_sphere(m, R))(H_state).simplify_full())
