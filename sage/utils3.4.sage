def L_polar(m, V):
    def Lagrangian(local):
        r, phi = coordinate(local)
        rdot, phidot = velocity(local)
        T = 1 / 2 * m * (square(rdot) + square(r * phidot))
        return T - V(r)

    return Lagrangian
