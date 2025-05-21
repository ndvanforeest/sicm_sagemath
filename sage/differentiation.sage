load(
    "functions.sage",
    "tuples.sage",
)


@Func
def D(expr):
    return lambda t: diff(expr(t), t)
    # return derivative(expr, t)


def Jacobian(F):
    def f(args, vrs):
        if isinstance(args, (list, tuple)):
            args = vector(args)
        if isinstance(vrs, (list, tuple)):
            vrs = vector(vrs)
        subs = {
            v: var(f"v{id(v)}", domain=RR)
            for v in args.list()
            if not v.is_symbol()
        }
        result = jacobian(F(args.subs(subs)), vrs.subs(subs).list())
        inverse_subs = {v: k for k, v in subs.items()}
        return result.subs(inverse_subs)

    return f


def gradient(F):
    return lambda v: Jacobian(F)(v, v).T


def Hessian(F):
    return lambda v: compose(gradient, gradient)(F)(v)


@Func
def partial(f, slot):
    def wrapper(local):
        if slot == 0:
            selection = [time(local)]
        elif slot == 1:
            selection = coordinate(local)
        elif slot == 2:
            selection = velocity(local)
        return Jacobian(f)(local, selection)

    return wrapper
