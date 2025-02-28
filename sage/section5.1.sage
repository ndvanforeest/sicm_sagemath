load(
    "utils.sage",
    "utils1.6.sage",
    "utils3.1.sage",
)

var("t x y", domain="real")

var("r, phi, p_r, p_phi", domain="real")
assume(r > 0)

q = vector([r, phi])
p = vector([p_r, p_phi])
local = up(t, q, p)

show((F_to_CH(p_to_r))(local)[0].simplify_full())
show((F_to_CH(p_to_r))(local)[1].simplify_full())
show((F_to_CH(p_to_r))(local)[2].simplify_full())

def H_central(m, V):
    def f(state):
        x, p = coordinate(state), momentum(state)
        return square(p) / (2 * m) + V(sqrt(square(x)))

    return f

var("r, phi, p_r, p_phi", domain="real")
assume(r > 0)
var("m", domain="positive")

show(
    Compose(H_central(m, function("V")), (F_to_CH(p_to_r)))(local)
    .simplify_full()
    .expand()
)

var("q_x q_y v_x v_y p_x p_y", domain="real")
q = vector([q_x, q_y])
v = vector([v_x, v_y])
p = vector([p_x, p_y])
local = up(t, q, p)

show(F_to_K(translating(v))(local))

def H_free(m):
    def f(state):
        return square(momentum(state)) / (2 * m)

    return f


def H_prime():
    return Sum(
        Compose(H_free(m), F_to_CH(translating(v))), F_to_K(translating(v))
    )


show(H_prime()(local))
