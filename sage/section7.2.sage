load("utils6.4.sage")

var("t", domain="real")

def H0(alpha):
    def f(state):
        p = momentum(state)[0, 0]
        return square(p) / 2 / alpha

    return f


def H1(beta):
    def f(state):
        theta = coordinate(state)[0, 0]
        return -beta * cos(theta)

    return f

def H_pendulum_series(alpha, beta, epsilon):
    def f(state):
        return H0(alpha)(state) + epsilon * H1(beta)(state)

    return f

def W(alpha, beta):
    def f(state):
        theta = coordinate(state)[0, 0]
        p = momentum(state)[0, 0]
        return -alpha * beta * sin(theta) / p

    return f

_ = var("theta p", domain="real")
_ = var("alpha beta", domain="positive")
epsilon = var("epsilon", domain="positive")


state = up(t, column_matrix([theta]), row_matrix([p]))
show(Lie_derivative(W(alpha, beta))(H0(alpha))(state))
show(H1(beta)(state))

show(
    Lie_transform(W(alpha, beta), epsilon)(
        H_pendulum_series(alpha, beta, epsilon)
    )(state, 2).expand()
)

show(Lie_transform(W(alpha, beta), epsilon)(time)(state, 2))

def solution0(alpha, beta):
    def f(t):
        def g(state0):
            t0 = time(state0)
            theta0 = coordinate(state0)[0, 0]
            p0 = momentum(state0)[0, 0]

            return up(
                t,
                column_matrix([theta0 + p0 * (t - t0) / alpha]),
                row_matrix([p0]),
            )

        return g

    return f

_ = var("delta_t", domain=RR, latex_name=r"\d t")
sol0 = solution0(alpha, beta)(delta_t)(state)
show(sol0)

def C(alpha, beta, epsilon, order):
    def f(state):
        return up(
            t,
            column_matrix(
                Lie_transform(W(alpha, beta), epsilon)(coordinate)(
                    state, order
                ).list()
            ),
            row_matrix(
                Lie_transform(W(alpha, beta), epsilon)(momentum)(
                    state, order
                ).list()
            ),
        )

    return f

show(C(alpha, beta, epsilon, 2)(state))

def C_inv(alpha, beta, epsilon, order):
    return C(alpha, beta, -epsilon, order)

def solution(epsilon, order):
    def f(alpha, beta):
        def g(delta_t):
            return compose(
                C(alpha, beta, epsilon, order),
                solution0(alpha, beta)(delta_t),
                C_inv(alpha, beta, epsilon, order),
            )

        return g

    return f

order = 1
sol = solution(epsilon, order)(alpha, beta)(delta_t)(state)

show(coordinate(sol)[0,0].expand())
