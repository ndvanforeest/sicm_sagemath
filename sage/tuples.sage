"""
This is a copy of tuples.py from
https://github.com/jtauber/functional-differential-geometry.
"""

from sage.structure.element import Matrix, Vector


class Tuple:
    def __init__(self, *components):
        self._components = components

    def __getitem__(self, index):
        return self._components[index]

    def __len__(self):
        return len(self._components)

    def __eq__(self, other):
        if (
            isinstance(other, self.__class)
            and self._components == other._components
        ):
            return True
        else:
            return False

    def __ne__(self, other):
        return not (self.__eq__(other))

    def __add__(self, other):
        if isinstance(self, Tuple):
            if not isinstance(other, self.__class__) or len(self) != len(
                other
            ):
                raise TypeError("can't add incompatible Tuples")
            else:
                return self.__class__(
                    *(
                        s + o
                        for (s, o) in zip(self._components, other._components)
                    )
                )
        else:
            return self + other

    def __iadd__(self, other):
        return self + other

    def __neg__(self):
        return self.__class__(*(-s for s in self._components))

    def __sub__(self, other):
        return self + (-other)

    def __isub__(self, other):
        return self - other

    def __call__(self, **kwargs):
        return self.__class__(
            *(
                (c(**kwargs) if isinstance(c, Expr) else c)
                for c in self._components
            )
        )

    def subs(self, args):
        # substitute variables with args
        return self.__class__(*(c.subs(args) for c in self._components))

    def list(self):
        "convert tuple and its components to one list."
        result = []
        for comp in self._components:
            if isinstance(comp, (Tuple, Matrix, Vector)):
                result.extend(comp.list())
            else:
                result.append(comp)
        return result

    def derivative(self, var):
        "Compute the derivative of all components and put the result in a tuple."
        return self.__class__(
            *[derivative(comp, var) for comp in self._components]
        )

class UpTuple(Tuple):
    def __repr__(self):
        return "up({})".format(", ".join(str(c) for c in self._components))

    def _latex_(self):
        "Print up tuples vertically."
        res = r"\begin{array}{c}"
        for comp in self._components:
            res += r"\begin{array}{c}"
            res += latex(comp)
            res += r"\end{array}"
            res += r" \\"
        res += r"\end{array}"
        return res

class DownTuple(Tuple):
    def __repr__(self):
        return "down({})".format(", ".join(str(c) for c in self._components))

    def _latex_(self):
        "Print down tuples horizontally."
        res = r"\begin{array}{c}"
        for comp in self._components:
            res += r"\begin{array}{c}"
            res += latex(comp)
            res += r"\end{array}"
            res += r" & "
        res += r"\end{array}"
        return res

up = UpTuple
down = DownTuple

up._dual = down
down._dual = up

def ref(tup, *indices):
    if indices:
        return ref(tup[indices[0]], *indices[1:])
    else:
        return tup


def component(*indices):
    def _(tup):
        return ref(tup, *indices)

    return _
