load("utils.sage")

q = vector(
    [
        literal_function("q_x"),
        literal_function("q_y"),
    ]
)

show(q(t=t))

qdot = derivative(q, t)
show(qdot)

up_tup = qv_to_state(q, qdot)
show(velocity(up_tup))

show(make_named_coordinates(["\\phi", r"\theta"]))
show(make_named_velocities(["\\phi", r"\theta"]))

show(make_velocities("q", dim=2))

space = make_named_space(["\\phi", "\\theta"])
show(space)

show(Gamma(q)(t))
