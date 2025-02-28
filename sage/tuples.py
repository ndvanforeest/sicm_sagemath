"""
This is a copy of tuples.py from https://github.com/jtauber/functional-differential-geometry
See ``tuples.rst`` in that repo for an explanation.
"""

# from symbolic import Expr, Add, Sub, Mul


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
            # if type(self) == type(other) and self._components == other._components:
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


class UpTuple(Tuple):
    def __repr__(self):
        return "up({})".format(", ".join(str(c) for c in self._components))


class DownTuple(Tuple):
    def __repr__(self):
        return "down({})".format(", ".join(str(c) for c in self._components))


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
