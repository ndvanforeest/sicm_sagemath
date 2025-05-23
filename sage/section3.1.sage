load("utils3.1.sage")

t = var("t", domain="real")

q = column_path([literal_function("q_x"), literal_function("q_y")])
p = row_path([literal_function("p_x"), literal_function("p_y")])

show(p(t))

H_state = qp_to_H_state_path(q, p)(t)
show(H_state)

V = Function(lambda x: function("V")(*x.list()))

H = H_rectangular
show(H(m, V)(H_state))

show(partial(H(m, V), 1)(H_state))

show(Hamiltonian_to_state_derivative(H(m, V))(H_state))

show(Hamilton_equations(H(m, V))(q, p)(t))

res = Lagrangian_to_Hamiltonian(L_central_rectangular(m, V))(H_state)
show(res)

show(res.simplify_full())

var("m g l")
q = column_path([literal_function("theta")])
p = row_path([literal_function("p")])

# space = make_named_space(["\\theta"])
H_state = qp_to_H_state_path(q, p)(t)
show(Lagrangian_to_Hamiltonian(L_planar_pendulum(m, g, l))(H_state))

q = column_path([literal_function("q_x"), literal_function("q_y")])
p = row_path([literal_function("p_x"), literal_function("p_y")])
H_state = qp_to_H_state_path(q, p)(t)
show(Lagrangian_to_Hamiltonian(L_Henon_Heiles(m))(H_state))

def L_sphere(m, R):
    def Lagrangian(local):
        theta, phi = coordinate(local).list()
        thetadot, phidot = velocity(local).list()
        return 1 / 2 * m * R ^ 2 * (
            square(thetadot) + square(phidot * sin(theta))
        )

    return Lagrangian


var("R", domain="positive")

q = column_path([literal_function("theta"), literal_function("phi")])
p = row_path([literal_function("p_x"), literal_function("p_y")])
H_state = qp_to_H_state_path(q, p)(t)
show(Lagrangian_to_Hamiltonian(L_sphere(m, R))(H_state).simplify_full())
