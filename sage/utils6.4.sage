load(
    "utils.sage",
    "utils1.6.sage",
    "utils3.1.sage",
    "utils3.2.sage",
    "utils3.4.sage",
)

var("t x y", domain="real")

def Lie_derivative(H):
    def f(F):
        return Poisson_bracket(F, H)

    return f

def Lie_transform(H, t, n):
    lie = Lie_derivative(H)

    def outer(func):
        def inner(local):
            term = func
            res = term(local)
            for i in range(1, n + 1):
                term = lie(term)
                res += t ^ i * term(local) / factorial(i)
            return res

        return inner

    return outer
