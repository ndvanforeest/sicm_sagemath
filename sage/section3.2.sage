load("utils3.2.sage")

t = var("t", domain="real")

q = column_matrix([var("q_x"), var("q_y")])
p = row_matrix([var("p_x"), var("p_y")])
sigma = up(t, q, p)
H = state_function("H")

show(Poisson_bracket(coordinate, H)(sigma))
show(Poisson_bracket(momentum, H)(sigma))

show(Poisson_bracket(H, H)(sigma))

F = state_function("F")
G = state_function("G")

show((Poisson_bracket(F, G) + Poisson_bracket(G, F))(sigma))

show(
    (
        Poisson_bracket(F, G + H)
        - Poisson_bracket(F, G)
        - Poisson_bracket(F, H)
    )(sigma)
)

constant = Function(lambda H_state: function("c")())

show(Jacobian(constant)(sigma, sigma))

show(
    (Poisson_bracket(F, constant * G) - constant * Poisson_bracket(F, G))(
        sigma
    ).simplify_full()
)

jacobi = (
    Poisson_bracket(F, Poisson_bracket(G, H))
    + Poisson_bracket(G, Poisson_bracket(H, F))
    + Poisson_bracket(H, Poisson_bracket(F, G))
)

show(jacobi(sigma).simplify_full())

def f(H_state):
    return function("f")(
        *coordinate(H_state).list(), *momentum(H_state).list()
    )

show(diff(f(sigma), time(sigma)))

V = Function(lambda q: function("V")(*q.list()))

var(m, domain="positive")

H = H_rectangular(m, V)

q = column_matrix([var("q")])
p = row_matrix([var("p")])
sigma = up(t, q, p)

show(Poisson_bracket(f, H)(sigma).expand())
