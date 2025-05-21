load("utils1.6.sage")

def f_bar(q_prime):
    q = lambda t: compose(F, Gamma(q_prime))(t)
    return lambda t: Gamma(q)(t)

def osculating_path(local):
    t = time(local)
    q = coordinate(local)

    def wrapper(t_prime):
        res = q
        dt = 1
        factorial = 1
        for k in range(2, len(local)):
            factorial *= k
            dt *= t_prime - t
            res += local[k] * dt / factorial
        return res

    return wrapper

def Gamma_bar(f_bar):
    def wrapped(local):
        t = time(local)
        q_prime = osculating_path(local)
        return f_bar(q_prime)(t)

    return wrapped

def F_to_C(F):
    def C(local):
        n = len(local)

        def f_bar(q_prime):
            q = lambda t: compose(F, Gamma(q_prime))(t)
            return lambda t: Gamma(q, n)(t)

        return Gamma_bar(f_bar)(local)

    return C

@Func
def Dt(F):
    def DtF(local):
        n = len(local)

        def DF_on_path(q):
            return D(lambda t: F(Gamma(q, n - 1)(t)))

        return Gamma_bar(DF_on_path)(local)

    return lambda state: DtF(local)

def Euler_Lagrange_operator(L):
    return lambda local: (Dt(partial(L, 2)) - partial(L, 1))(local)
