import numpy as np

load("functions.sage", "differentiation.sage", "tuples.sage")

def L_free_particle(mass):
    def Lagrangian(local):
        v = velocity(local)
        return 1 / 2 * mass * square(v)

    return Lagrangian

@Func
def literal_function(name):
    return lambda t: function(name)(t)

@Func
def column_path(lst):
    #return lambda t: vector([l(t) for l in lst])
    return lambda t: column_matrix([l(t) for l in lst])

def Gamma(q, n=3):
    # if isinstance(q, np.ndarray):
    #     q = vector(q.tolist()) # todo, is this still needed?

    if n < 2:
        raise ValueError("n must be > 1")
    Dq = [q]
    for k in range(2, n):
        Dq.append(D(Dq[-1]))
    return lambda t: up(t, *[v(t) for v in Dq])

time = Function(lambda local: local[0])
coordinate = Function(lambda local: local[1])
velocity = Function(lambda local: local[2])

from sage.symbolic.integration.integral import definite_integral

def integral_latex_format(*args):
    expr, var, a, b = args
    return (
        fr"\int_{{{a}}}^{{{b}}} "
        + latex(expr)
        + r"\, \textrm{d}\,"
        + latex(var)
    )


definite_integral._print_latex_ = integral_latex_format

def Lagrangian_polynomial(ts, qs):
    return RR['x'].lagrange_polynomial(list(zip(ts, qs)))

def L_harmonic(m, k):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        return (1 / 2) * m * square(v) - (1 / 2) * k * square(q)

    return Lagrangian
