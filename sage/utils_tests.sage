load("utils.sage")

q = vector(
    [
        literal_function("q_x"),
        literal_function("q_y"),
        literal_function("q_z"),
    ]
)

show(q(t=t))

show(q.diff(t))

show(make_named_coordinates(["q", "r"]))
show(make_named_velocities(["q", "r"]))
show(make_named_coordinates(["\\phi", r"\theta"]))
show(make_named_velocities(["\\phi", r"\theta"]))

show(make_coordinates("q", dim=3))
show(make_velocities("q", dim=3))

show(make_space("q", dim=2))

space = make_named_space(["\\phi", "\\theta"])
show(space)

var("a b c x y", domain="real")

F = a * square(x) + b * x + c
show(diff(F, x))
show(diff(F, x)(x=0))
show(diff(F, x, 2))
show(diff(F, x, x))

M = matrix([[3, 4], [4, 5]])
b = vector([8, 9])
v = vector([x, y])
F = 1 / 2 * v * M * v + b * v + c
show(F.args())

show(jacobian(F, (x, y)))
show(jacobian(F, v))
show(jacobian(F, (x, y))(x=0, y=0))
show(jacobian(F(x=x, y=y), (x, y)))
show(jacobian(F(x=x, y=y), (x, y))(x=0, y=0))

show(jacobian(jacobian(F, (x, y)), (x, y)))

def F(v):
    return 1 / 2 * v * M * v + b * v + c

show(jacobian(F(v), (x, y)))
show(jacobian(F(v), v))
# show(jacobian(F, v))  # This does not work

U = function("U")
show(jacobian(U(*v), (x, y)))

def U(q):
    return function("U")(*q)

show(jacobian(U(v), (x, y)))

def gradient(F, v):
    cds = make_coordinates(f"qq_{id(F)}", dim=len(v))
    if isinstance(F, sage.symbolic.function_factory.SymbolicFunction):
        deriv = jacobian(F(*cds), cds)  # Unpack coordinates if F is symbolic
    else:
        deriv = jacobian(F(cds), cds)  # Otherwise, call F with a vector
    return deriv.subs(dict(zip(cds, v)))

show(gradient(F, v))

show(gradient(F, v)(x=0, y=0))
show(gradient(F, v).subs({v[0]: 0, v[1]: 0}))
show(gradient(F, [0, 0]))

q = vector([literal_function("q_1"), literal_function("q_2")])
show(gradient(F, q))

def U(q):
    return function("U")(*q)

show(gradient(U, q))

show(hessian(F, q))

def U(q):
    return function("U")(*q)

show(hessian(U, q))

def L_harmonic(m, k):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        return (1 / 2) * m * square(v) - (1 / 2) * k * square(q)

    return Lagrangian


var('k m', domain="positive")
L = L_harmonic(m, k)

space = make_space("x", dim=2)
show(space)
show(partial(L, 1)(space))
show(partial(L, 2)(space))

def L_generic(m, U):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        return (1 / 2) * m * square(v) - U(q)

    return Lagrangian

def U(q):
    return function("U")(*q)


L_gen = L_generic(m, U)
show(partial(L_gen, 1)(space))

q = vector([literal_function("q_1"), literal_function("q_2")])
show(partial(L, 1)(Gamma(q)(t)))
show(partial(L_gen, 1)(Gamma(q)(t)))

show(D(q))

f = function('f')
g = function('g')
h = function('h')

show(Sum(f, g, h)(x))
show(Product(f, g, h)(x))
show(Compose(f, g, h)(x))
show(diff(Sum(f, g, h)(x), x))
show(diff(Product(f, g, h)(x), x))
show(diff(Compose(f, g, h)(x), x))

f = lambda x: sin(x)
g = lambda x: x ^ 2 + 1
h = lambda x: 3 * x


show(Sum(f, g, h)(x))
show(Product(f, g, h)(x))
show(Compose(f, g, h)(x))
show(diff(Sum(f, g, h)(x), x))
show(diff(Product(f, g, h)(x), x))
show(diff(Compose(f, g, h)(x), x))

show(Gamma(q)(t))

show(velocity(Gamma(q)(t)))
