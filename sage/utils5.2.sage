def vector_to_up(vec):
    n = len(vec)
    dim = n // 2
    return up(vec[0], vec[1 : dim + 1], vec[dim + 1 : n + 1])


def up_to_vector(tup):
    return vector([time(tup), *coordinate(tup), *velocity(tup)])

def to_vector_map(C):
    "Lift the structure map C such that it maps a vector to a vector."
    return Compose(up_to_vector, C, vector_to_up)

def D_as_matrix(C):
    def f(state):
        vec = up_to_vector(state)
        vmap = to_vector_map(C)
        return jacobian(vmap(vec), vec).simplify_full()

    return f

def symplectic_unit(n):
    I = identity_matrix(n)
    return block_matrix([[zero_matrix(n), I], [-I, zero_matrix(n)]])

def is_symplectic_matrix(M):
    n = M.nrows()
    J = symplectic_unit(n // 2)
    M = M.expand()
    res = (M * J * M.transpose()).simplify_full()
    if res == J:
        return True
    print(res - J)
    return False

def qp_submatrix(M):
    return M[1:, 1:]

def is_symplectic_transform(C):
    return lambda state: Compose(
        is_symplectic_matrix, qp_submatrix, D_as_matrix(C)
    )(state)
