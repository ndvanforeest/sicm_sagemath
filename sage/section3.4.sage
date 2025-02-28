import numpy as np

load(
    "utils.sage",
    "utils3.1.sage",
)

var("t x y", domain="real")

var("m")

V = function("V")

def L_polar(m, V):
    def Lagrangian(local):
        r, phi = coordinate(local)
        rdot, phidot = velocity(local)
        T = 1 / 2 * m * (square(rdot) + square(r * phidot))
        return T - V(r)

    return Lagrangian

space = make_named_space(["r", "\\phi"])
show(space)

q = vector([literal_function("r"), literal_function("\\phi")])
p = vector([literal_function(r"p_r"), literal_function(r"p_{\phi}")])
H_state = qp_to_H_state_path(q, p)(t)
show(H_state)

H = Lagrangian_to_Hamiltonian(L_polar(m, V))  # (H_state)
show(H(H_state))

show(H(H_state).expand())

HE = Hamilton_equations(Lagrangian_to_Hamiltonian(L_polar(m, V)))(q, p)
show(HE)
