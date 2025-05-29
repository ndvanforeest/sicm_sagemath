load("tuples.sage")

class Function:
    def __init__(self, func):
        self._func = func

    def __call__(self, *args):
        return self._func(*args)

    def __add__(self, other):
        return Function(lambda *args: self(*args) + other(*args))

    def __neg__(self):
        return Function(lambda *args: -self(*args))

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, Function):
            return Function(lambda *args: self(*args) * other(*args))
        return Function(lambda *args: other * self(*args))

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exponent):
        if exponent == 0:
            return Function(lambda x: 1)
        else:
            return self * (self ** (exponent - 1))

def Func(f):
    def wrapper(*args, **kwargs):
        return Function(f(*args, **kwargs))

    return wrapper

@Func
def compose(*funcs):
    if len(funcs) == 1:
        return lambda x: funcs[0](x)
    return lambda x: funcs[0](compose(*funcs[1:])(x))

identity = Function(lambda x: x)

sin = Function(lambda x: sage.functions.trig.sin(x))
cos = Function(lambda x: sage.functions.trig.cos(x))

from functools import singledispatch


@singledispatch
def _square(x):
    raise TypeError(f"Unsupported type: {type(x)}")


@_square.register(int)
@_square.register(float)
@_square.register(Expression)
@_square.register(Integer)
def _(x):
    return x ^ 2


@_square.register(Vector)
@_square.register(list)
@_square.register(tuple)
def _(x):
    v = vector(x)
    return v.dot_product(v)


@_square.register(Matrix)
def _(x):
    if x.ncols() == 1:
        return (x.T * x)[0, 0]
    elif x.nrows() == 1:
        return (x * x.T)[0, 0]
    else:
        raise TypeError(
            f"Matrix must be a row or column vector, got shape {x.nrows()}×{x.ncols()}"
        )


square = Function(lambda x: _square(x))

function = sage.symbolic.function_factory.function

cube = Function(lambda x: x * square(x))
