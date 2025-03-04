import numpy as np
from tuples import up, Tuple  # See below why this import.

var('t', domain="real")


def literal_function(name):
    return function(name, nargs=1, print_latex_func=print_lit_f_to_latex)(t)


def print_lit_f_to_latex(name, *args):
    return name


def qv_to_state(q, v):
    return up(t, q, v)


def qv_to_state_path(q, v):
    def f(t):
        return up(t, q(t=t), v(t=t))

    return f


def make_named_coordinates(coordinate_names, latex_names=None):
    vars = []
    if latex_names == None:
        latex_names = [name for name in coordinate_names]

    stripped = [f'{name.lstrip(r"\\")}' for name in coordinate_names]
    for name, latex in zip(stripped, latex_names):
        q = var(name, latex_name=f"{latex}", domain='real')
        vars.append(q)
    # return column_matrix(vars)
    return vector(vars)


def make_named_velocities(coordinate_names, latex_names=None):
    names = [f"{name}dot" for name in coordinate_names]
    if latex_names == None:
        latex_names = [fr"\dot {name}" for name in coordinate_names]
    return make_named_coordinates(names, latex_names)


def make_coordinates(coordinate_name, dim):
    names = [f"{coordinate_name}_{i}" for i in range(1, dim + 1)]
    return make_named_coordinates(names)


def make_velocities(coordinate_name, dim):
    names = [f"{coordinate_name}_{i}" for i in range(1, dim + 1)]
    return make_named_velocities(names)


def make_named_space(coordinate_names):
    coordinates = make_named_coordinates(coordinate_names)
    velocities = make_named_velocities(coordinate_names)
    return qv_to_state(coordinates, velocities)


def make_space(coordinate_name, dim):
    coordinates = make_coordinates(coordinate_name, dim)
    velocities = make_velocities(coordinate_name, dim)
    return qv_to_state(coordinates, velocities)


def square(x):
    return x * x


def gradient(F, v):
    cds = make_coordinates(f"q_{id(F)}", dim=len(v))
    deriv = jacobian(F(cds), cds)
    return vector(deriv.subs(dict(zip(cds, v))))


def hessian(F, v):
    cds = make_coordinates(f"q_{id(F)}", dim=len(v))
    hes = jacobian(jacobian(F(cds), cds), cds)
    return matrix(hes.subs(dict(zip(cds, v))))


def partial(f, slot):
    def wrapper(local):
        space = make_space(f"q_{id(f)}_{slot}", dim=len(coordinate(local)))
        if slot == 0:
            selection = [time(space)]
        elif slot == 1:
            selection = coordinate(space)
        elif slot == 2:
            selection = velocity(space)
        deriv = jacobian(f(space), selection)
        return deriv.subs(
            {
                t: time(local),
                **dict(zip(coordinate(space), coordinate(local))),
                **dict(zip(velocity(space), velocity(local))),
            }
        )

    return wrapper


def D(expr):
    "Derivative wrt time t."
    if isinstance(expr, Tuple):
        return up(
            D(time(expr)),
            D(coordinate(expr)),
            D(velocity(expr)),
        )
    return derivative(expr, t)


def Sum(*funcs):
    if len(funcs) == 1:
        return lambda x: funcs[0](x)
    return lambda x: funcs[0](x) + Sum(*funcs[1:])(x)


def Product(*funcs):
    if len(funcs) == 1:
        return lambda x: funcs[0](x)
    return lambda x: funcs[0](x) * Product(*funcs[1:])(x)


def Compose(*funcs):
    if len(funcs) == 1:
        return lambda x: funcs[0](x)
    return lambda x: funcs[0](Compose(*funcs[1:])(x))


def Min(func):
    return lambda x: -func(x)


def Gamma(q):
    q = vector(q)
    v = D(q)
    return qv_to_state_path(q, v)


def time(local):
    return local[0]


def coordinate(local):
    return local[1]


def velocity(local):
    return local[2]


def rotation_matrix(axis, theta):
    """
    Return the 3x3 rotation matrix for a rotation of angle theta (in radians)
    about the given axis. The axis is specified as an iterable of 3 numbers.
    """
    # Convert the axis to a normalized vector
    axis = vector(axis).normalized()
    x, y, z = axis
    c = cos(theta)
    s = sin(theta)
    t = 1 - c  # common factor

    # Construct the rotation matrix using Rodrigues' formula
    R = matrix(
        [
            [c + x**2 * t, x * y * t - z * s, x * z * t + y * s],
            [y * x * t + z * s, c + y**2 * t, y * z * t - x * s],
            [z * x * t - y * s, z * y * t + x * s, c + z**2 * t],
        ]
    )
    return R
