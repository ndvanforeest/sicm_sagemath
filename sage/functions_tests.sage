load("functions.sage")

def f(x):
    return 5 * x


F = Function(lambda x: f(x))

V = Function(lambda x: function("V")(x))

x, y = var("x y", domain = RR)

show((square)(x + y).expand())

show((square + square)(x + y))

show((square * square)(x))

show((sin + cos)(x))

show((square + V)(x))

hh = compose(square, sin)
show((hh + hh)(x))

show((2 * (sin * cos)(x) - sin(2 * x)).simplify_full())

show(diff(-compose(square, cos)(x), x))
show(integrate((2 * sin * cos)(x), x))

U = Function(lambda x: function("U")(x))
V = Function(lambda x: function("V")(x))

show((U + V)(x))
show((V + V)(x))
show((V(U(x))))
show((compose(V, U)(x)))

def f(x):
    def g(y):
        return x * y ^ 2

    return g

show(f(3)(5))

@Func
def f(x):
    def g(y):
        return x * y ^ 2

    return g

show((f(3) + f(2))(4))

def f(x):
    def g(y):
        return x * y ^ 2

    return Function(lambda y: g(y))

show((f(3) + f(2))(4))
