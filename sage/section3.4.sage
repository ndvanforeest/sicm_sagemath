load("utils3.1.sage")

t = var("t", domain="real")



var("m")

V = function("V")


def L_polar(m, V):
    def Lagrangian(local):
        r, phi = coordinate(local).list()
        rdot, phidot = velocity(local).list()
        T = 1 / 2 * m * (square(rdot) + square(r * phidot))
        return T - V(r)

    return Lagrangian

q = column_path([literal_function("r"), literal_function("phi")])
p = row_path([literal_function("p_r"), literal_function("p_phi")])
H_state = qp_to_H_state_path(q, p)(t)
show(H_state)

H = Lagrangian_to_Hamiltonian(L_polar(m, V))
show(H(H_state).simplify_full())

HE = Hamilton_equations(Lagrangian_to_Hamiltonian(L_polar(m, V)))(q, p)(t)
show(HE)
