load("utils1.6.sage")


def Hamilton_equations(Hamiltonian):
    def f(q, p):
        state_path = qp_to_H_state_path(q, p)
        return D(state_path) - compose(
            Hamiltonian_to_state_derivative(Hamiltonian), state_path
        )

    return f


def qp_to_H_state_path(q, p):
    def f(t):
        return up(t, q(t), p(t))

    return f


def transpose(M):
    return M.T


def row_path(lst):
    return lambda t: transpose(column_path(lst)(t))


def row_matrix(lst):
    return transpose(column_matrix(lst))


def Hamiltonian_to_state_derivative(Hamiltonian):
    def f(H_state):
        return up(
            SR(1),
            partial(Hamiltonian, 2)(H_state).T,
            -partial(Hamiltonian, 1)(H_state),
        )

    return f


var("m")


def H_rectangular(m, V):
    def f(state):
        q, p = coordinate(state), momentum(state)
        return square(p) / 2 / m + V(q)

    return f


momentum = Function(lambda H_state: H_state[2])


def Legendre_transform(F):
    def G(w):
        zeros = {var(f"v_{i}"): 0 for i in range(w.ncols())}
        b = gradient(F)(list(zeros.keys())).subs(zeros)
        M = Hessian(F)(list(zeros.keys())).subs(zeros)
        v = M.solve_right(w.T - b)
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
