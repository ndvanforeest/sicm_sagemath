load("utils3.2.sage")

@Func
def Lie_derivative(H):
    def f(F):
        return Poisson_bracket(F, H)

    return f

def Lie_transform(H, t):
    lie = Lie_derivative(H)

    def outer(func):
        def inner(local, n):
            term = func
            factor = 1
            res = factor * term(local)
            for i in range(1, n + 1):
                term = lie(term)
                factor *= t / i
                res += factor * term(local)
            return res

        return inner

    return outer
