load(
    "utils1.6.sage",
    "utils3.1.sage",
)

def F_to_CH(F):
    M = partial(F, 1)

    def f(state):
        return up(time(state), F(state), M(state).solve_left(momentum(state)))

    return f

def F_to_K(F):
    M = partial(F, 1)

    def f(state):
        p = M(state).solve_left(momentum(state))
        return -p * partial(F, 0)(state)

    return f

def translating(v):
    def f(state):
        return coordinate(state) + v * time(state)

    return f
