load("utils5.1.sage")

t = var("t", domain="real")

var("r, phi, p_r, p_phi", domain="real")
assume(r > 0)

q = column_matrix([r, phi])
p = row_matrix([p_r, p_phi])
state = up(t, q, p)

show((F_to_CH(p_to_r))(state)[0].simplify_full())
show((F_to_CH(p_to_r))(state)[1].simplify_full())
show((F_to_CH(p_to_r))(state)[2].simplify_full())

def H_central(m, V):
    def f(state):
        x, p = coordinate(state), momentum(state)
        return square(p) / (2 * m) + V(sqrt(square(x)))

    return f

var("m", domain="positive")



show(
    compose(H_central(m, function("V")), (F_to_CH(p_to_r)))(state)
    .simplify_full()
    .expand()
)

var("q_x q_y v_x v_y p_x p_y", domain="real")
q = column_matrix([q_x, q_y])
v = column_matrix([v_x, v_y])
p = row_matrix([p_x, p_y])
state = up(t, q, p)

show(F_to_K(translating(v))(state)[0, 0])

def H_free(m):
    def f(state):
        return square(momentum(state)) / (2 * m)

    return f


def H_prime():
    return compose(H_free(m), F_to_CH(translating(v))) + F_to_K(translating(v))



show(H_prime()(state)[0, 0])
