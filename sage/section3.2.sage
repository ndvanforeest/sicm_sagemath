load(
    "utils.sage",
    "utils1.6.sage",
    "utils3.1.sage",
    "utils3.2.sage",
)

var("t x y", domain="real")

def Poisson_bracket(F, G):
    def f(state):
        left = partial(F, 1)(state) * partial(G, 2)(state).T
        right = partial(F, 2)(state) * partial(G, 1)(state).T
        return (left - right).simplify_full()

    return f

space = make_space("q", dim=2)

def F(local):
    return function("f")(time(local), *coordinate(local), *velocity(local))

show(Poisson_bracket(F, F)(space))

def G(local):
    return function("g")(time(local), *coordinate(local), *velocity(local))

show(Sum(Poisson_bracket(F, G), Poisson_bracket(G, F))(space))

def H(local):
    return function("h")(time(local), *coordinate(local), *velocity(local))

show(Poisson_bracket(F, Sum(G,H))(space) == Sum(Poisson_bracket(F, G), Poisson_bracket(F, H))(space))

def constant(local):
    return function("c")()

show(diff(constant(space), time(space)))
show(diff(constant(space), *coordinate(space)))

show(
    Poisson_bracket(F, Product(constant, G))(space)
    == Product(constant, Poisson_bracket(F, G))(space)
)

jacobi = Compose(
    simplify,
    Sum(
        Poisson_bracket(F, Poisson_bracket(G, H)),
        Poisson_bracket(G, Poisson_bracket(H, F)),
        Poisson_bracket(H, Poisson_bracket(F, G)),
    ),
)
show(jacobi(space))

def F(local):
    return function("f")(*coordinate(local), *velocity(local))


show(diff(F(space), time(space)))

def V(q):
    return function("U")(*q)


var(m, domain="positive")

H = H_rectangular(m, V)

space = make_space("q", dim=1)
show(Poisson_bracket(F, H)(space).expand()[0, 0])
