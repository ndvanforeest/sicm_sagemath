load("utils6.4.sage")

var("t", domain="real")

f = function("f")
show(taylor(f(x), x, 0, 4))

f = sin
show(taylor(f(x), x, 0, 6))

# show(taylor(lambda x: sqrt(1 + x), x, 0, 6))  # this does not work
show(taylor((lambda x: sqrt(1 + x))(x), x, 0, 6))

n = var("n", domain="integer")
f = lambda x: (1 + x) ^ n
show(taylor(f(x), x, 0, 7))

_ = var("m k", domain="positive")


def H_harmonic(m, k):
    def f(state):
        return (
            square(momentum(state)) / (2 * m)
            + k * square(coordinate(state)) / 2
        )

    return f

H = H_harmonic(m, k)

x0, p0 = var("x0 p0", domain="real")
state = up(t, column_matrix([x0]), row_matrix([p0]))
lie = Lie_derivative(H)
show(lie(coordinate)(state))
show(lie(lie(coordinate))(state))

_ = var('dt', domain="real", latex_name=r"\d t")

show(Lie_transform(H, dt)(coordinate)(state, 4))
show(Lie_transform(H, dt)(momentum)(state, 4))
show(Lie_transform(H, dt)(H)(state, 4))

def V(q):
    return function("U")(q)

m = var("m", domain="positive")

def H_central_polar(m, V):
    def f(state):
        r, phi = coordinate(state).list()
        p_r, p_phi = momentum(state).list()
        T = 1 / 2 * square(p_r) / m + 1 / 2 * square(p_phi) / (m * square(r))
        return T + V(r)

    return f


H = H_central_polar(m, V)

_ = var("r phi p_r p_phi", domain="real")
assume(r > 0)
q = column_matrix([r, phi])
p = row_matrix([p_r, p_phi])
H_state = up(t, q, p)

show(H(H_state).expand())
show(partial(H, 1)(H_state))

res = Lie_transform(H, dt)(coordinate)(H_state, 3).expand()
show(res[0][0])
show(res[1][0])
