def transpose(M):
    return M.T


def simplify(expr):
    return expr.simplify_full()


def Poisson_bracket(F, G):
    def f(state):
        return Compose(
            simplify,
            Sum(
                Product(partial(F, 1), Compose(transpose, partial(G, 2))),
                Min(Product(partial(F, 2), Compose(transpose, partial(G, 1)))),
            ),
        )(state)

    return f
