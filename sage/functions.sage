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

def _square(x):
    if isinstance(
        x,
        (
            int,
            float,
            sage.symbolic.expression.Expression,
            sage.rings.integer.Integer,
        ),
    ):
        return x ^ 2
    elif isinstance(x, (Vector, list, tuple)):
        v = vector(x)
        return v.dot_product(v)
    elif isinstance(x, Matrix) and x.ncols() == 1:
        return (x.transpose() * x)[0, 0]
    else:
        raise TypeError(f"Unsupported type: {type(x)}")


square = Function(lambda x: _square(x))

sin = Function(lambda x: sage.functions.trig.sin(x))
cos = Function(lambda x: sage.functions.trig.cos(x))

function = sage.symbolic.function_factory.function
