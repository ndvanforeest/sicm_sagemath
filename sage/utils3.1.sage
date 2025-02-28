load(
    "utils.sage",
    "utils1.6.sage",
)

def Hamilton_equations(Hamiltonian):
    def f(q, p):
        state_path = qp_to_H_state_path(q, p)(t)
        return D(state_path) - Hamiltonian_to_state_derivative(Hamiltonian)(
            state_path
        )

    return f

def Hamiltonian_to_state_derivative(Hamiltonian):
    def f(H_state):
        return up(
            SR(1),
            vector(partial(Hamiltonian, 2)(H_state)),
            -vector(partial(Hamiltonian, 1)(H_state)),
        )

    return f

def qp_to_H_state_path(q, p):
    def f(t):
        return up(t, q(t=t), p(t=t))

    return f

def momentum(H_state):
    return H_state[2]

def H_rectangular(m, V):
    def f(state):
        q, p = coordinate(state), velocity(state)
        return square(p) / 2 / m + V(q)

    return f

def Legendre_transform(F):
    def G(w):
        zero = [0] * len(w)
        b = gradient(F, zero)
        M = hessian(F, zero)
        v = M.solve_left(w - b)
        return w * v - F(v)

    return G

def Lagrangian_to_Hamiltonian(Lagrangian):
    def f(H_state):
        t = time(H_state)
        q = coordinate(H_state)
        p = momentum(H_state)

        def L(qdot):
            return Lagrangian(up(t, q, qdot))

        return Legendre_transform(L)(p)

    return f
