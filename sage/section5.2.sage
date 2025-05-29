load("utils5.2.sage")

t = var("t", domain="real")

show(symplectic_unit(2))

J = symplectic_unit(2)
show(is_symplectic_matrix(J))

var("r, phi, p_r, p_phi", domain="real")
assume(r > 0)

r_phi = up(t, column_matrix([r, phi]), row_matrix([p_r, p_phi]))

vec = up_to_vector(r_phi)
show(vec)

show(up_to_vector(F_to_CH(p_to_r)(r_phi)).simplify_full())

show(to_vector_map(F_to_CH(p_to_r))(r_phi).simplify_full())

show(is_symplectic_transform(F_to_CH(p_to_r))(r_phi))

def F(local):
    t, q = time(local), coordinate(local)
    return vector([function("f")(t, *q), function("g")(t, *q)])

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
