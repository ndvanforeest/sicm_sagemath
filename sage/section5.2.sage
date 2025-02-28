load(
    "utils.sage",
    "utils1.6.sage",
    "utils3.1.sage",
    "utils5.1.sage",
    "utils5.2.sage",
)

var("t x y", domain="real")

var("r, phi, p_r, p_phi", domain="real")
assume(r > 0)

r_phi = up(t, vector([r, phi]), vector([p_r, p_phi]))

vec = up_to_vector(r_phi)
show(vec)
show(up_to_vector(vector_to_up(vec)))
show(vector_to_up(up_to_vector(r_phi)))

show(up_to_vector(F_to_CH(p_to_r)(r_phi)).simplify_full())

CH = lambda t: up_to_vector(F_to_CH(p_to_r)(r_phi))

res = diff(CH(t), t)
show(res)

DC_H = jacobian(to_vector_map(F_to_CH(p_to_r))(vec), vec).simplify_full()
show(DC_H)

show(symplectic_unit(2))

J = symplectic_unit(2)
show(is_symplectic_matrix(J))

show(is_symplectic_transform(F_to_CH(p_to_r))(r_phi))

_ = var("q_x q_y v_x v_y p_x p_y", domain="real")
xy = up(t, vector([x, y]), vector([p_x, p_y]))


def F(local):
    t = time(local)
    q = coordinate(local)
    return vector([function("f")(t, *q), function("g")(t, *q)])


show(is_symplectic_transform(F_to_CH(F))(xy))

def polar_canonical(alpha):
    def f(state):
        t = time(state)
        theta = coordinate(state)[0]
        I = momentum(state)[0]
        x = sqrt(2 * I / alpha) * sin(theta)
        p = sqrt(2 * I * alpha) * cos(theta)
        return up(t, vector([x]), vector([p]))

    return f


_ = var("theta", domain="real")
_ = var("I alpha", domain="positive")
theta_I = up(t, vector([theta]), vector([I]))
show(is_symplectic_transform(polar_canonical(alpha))(theta_I))

def a_non_canonical_transform(state):
    t = time(state)
    theta = coordinate(state)[0]
    I = momentum(state)[0]
    x = I * sin(theta)
    p = I * cos(theta)
    return up(t, vector([x]), vector([p]))


show(is_symplectic_transform(a_non_canonical_transform)(theta_I))
