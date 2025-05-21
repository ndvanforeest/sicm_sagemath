load(
    "utils6.4.sage",
)

var("t x y", domain="real")

f = function("f")
show(taylor(f(x), x, 0, 4))

f = sin
show(taylor(f(x), x, 0, 6))

# show(taylor(lambda x: sqrt(1 + x), x, 0, 6))  # this does not work
show(taylor((lambda x: sqrt(1 + x))(x), x, 0, 6))

n = var("n", domain="integer")
f = lambda x: (1 + x) ^ n
show(taylor(f(x), x, 0, 4))

var("m k", domain="positive")


def H_harmonic(m, k):
    def f(state):
        return (
            square(momentum(state)) / (2 * m)
            + k * square(coordinate(state)) / 2
        )

    return f

H = H_harmonic(m, k)

x0, p0 = var("x0 p0", domain="real")
space = up(t, vector([x0]), vector([p0]))
lie = Lie_derivative(H)
show(lie(coordinate)(space))

show(lie(lie(coordinate))(space))
show(lie(lie(lie(coordinate)))(space))

# show((t * lie(coordinate))(space)) # gives an error
show(t * lie(coordinate)(space))

_ = var('dt', domain="real")

F = Compose(lambda x: x[0], coordinate)
show(Lie_transform(H, dt, 4)(F)(space))
G = Compose(lambda x: x[0], momentum)
show(Lie_transform(H, dt, 4)(G)(space))

F = Compose(column_matrix, coordinate)
show(Lie_transform(H, dt, 4)(H)(space))

def V(q):
    return function("U")(q)

m = var("m", domain="positive")

H = Lagrangian_to_Hamiltonian(L_polar(m, V))

def H_central_polar(m, V):
    def f(state):
        r, phi = coordinate(state)
        p_r, p_phi = momentum(state)
        T = 1 / 2 * square(p_r) / m + 1 / 2 * square(p_phi) / (m * square(r))
        return T + V(r)

    return f


H = H_central_polar(m, V)

_ = var("r phi p_r p_phi", domain="real")
assume(r > 0)
q = vector([r, phi])
p = vector([p_r, p_phi])
H_state = up(t, q, p)

show(H(H_state).expand())
show(partial(H, 1)(H_state))

res = Lie_transform(H, dt, 3)(Compose(column_matrix, coordinate))(H_state).expand()
show(res[0][0])
show(res[1][0])
