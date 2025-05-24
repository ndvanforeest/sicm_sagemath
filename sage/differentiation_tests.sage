load("differentiation.sage")

var("t", domain="real")

_ = var("a b c x y", domain=RR)
M = matrix([[a, b], [b, c]])
b = vector([a, b])
v = vector([x, y])
F = 1 / 2 * v * M * v + b * v + c

show(F)

show(F.expand())

show(diff(F, x))

show(diff(F, [x, y]))

show(jacobian(F, [x, y]))

show(jacobian(F, v.list()))  # convert the column matrix to a list

def F(v):
    return 1 / 2 * v * M * v + b * v + c

show(diff(F(v), x)) # add the arguments to F
show(jacobian(F(v), v.list()))

T = up(t, t ^ 2, t ^ 3, sin(3 * t))
show(diff(T, t))

q = Function(lambda t: function("q")(t))

v = var("v", domain=RR)


def F(v):
    r, t = v.list()
    return 5 * r ^ 3 + 3 * t ^ 2 * r


show(Jacobian(F)([v, t], [t]))
show(Jacobian(F)([v, t], [v, t]))

q = Function(lambda t: function("q")(t))
v = [q(t), t]
show(Jacobian(F)(v, v))

show(gradient(F)(v))

U = Function(lambda x: function("U")(square(x)))
show(gradient(U)(v))

show(Hessian(F)(v))
