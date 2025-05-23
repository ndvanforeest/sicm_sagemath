load("utils3.1.sage")

@Func
def Poisson_bracket(F, G):
    def f(state):
        left = (partial(F, 1) * compose(transpose, partial(G, 2)))(state)
        right = (partial(F, 2) * compose(transpose, partial(G, 1)))(state)
        return (left - right).simplify_full()

    return f

@Func
def state_function(name):
    return lambda H_state: function(name)(
        time(H_state), *coordinate(H_state).list(), *momentum(H_state).list()
    )
