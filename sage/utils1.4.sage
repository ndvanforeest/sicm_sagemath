def L_free_particle(mass):
    def Lagrangian(local):
        v = velocity(local)
        return 1 / 2 * mass * square(v)

    return Lagrangian

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

def Lagrangian_action(L, q, t1, t2):
    return definite_integral(Compose(L, Gamma(q))(t), t, t1, t2)

def Lagrangian_polynomial(ts, qs):
    return RR['x'].lagrange_polynomial(list(zip(ts, qs)))

def L_harmonic(m, k):
    def Lagrangian(local):
        q = coordinate(local)
        v = velocity(local)
        return (1 / 2) * m * square(v) - (1 / 2) * k * square(q)

    return Lagrangian
