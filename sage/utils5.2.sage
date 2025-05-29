load("utils5.1.sage")

def symplectic_unit(n):
    I = identity_matrix(n)
    return block_matrix([[zero_matrix(n), I], [-I, zero_matrix(n)]])

def is_symplectic_matrix(M):
    n = M.nrows()
    J = symplectic_unit(n // 2)
    if isinstance(M, sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense):
        M = M.expand()
        res = (M * J * M.transpose()).simplify_full()
    else:
        res = M * J * M.transpose()
    if res == J:
        return True
    print(res - J)
    return False

def is_symplectic_transform(C):
    return lambda state: compose(
        is_symplectic_matrix, qp_submatrix, D_as_matrix(C)
    )(state)

def qp_submatrix(M):
    return M[1:, 1:]

def D_as_matrix(C):
    def f(state):
        vec = up_to_vector(state)
        vmap = to_vector_map(C)
        #return jacobian(vmap(vec), vec).simplify_full()
        return jacobian(vmap(state), vec).simplify_full()

    return f

def up_to_vector(state):
    return vector(
        [time(state), *coordinate(state).list(), *momentum(state).list()]
    )

def to_vector_map(C):
    return compose(up_to_vector, C)
