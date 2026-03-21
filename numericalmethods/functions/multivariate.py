from __future__ import annotations

from .base import Function
from ..linalg.vector import Vector


class MultiVariableFunction:
    """A function of two or more scalar (or :class:`~numericalmethods.linalg.Vector`) arguments.

    Supports **partial application**: pass ``None`` for any argument to defer it,
    returning a new :class:`~numericalmethods.functions.base.Function` (for a single
    free argument) or a new :class:`MultiVariableFunction` (for multiple free arguments).

    Parameters
    ----------
    function:
        A callable that accepts all positional scalar arguments.
    """

    def __init__(self, function) -> None:
        self.function = function

    def __call__(self, *args):
        # Unpack any Vector arguments into their scalar components
        unwrapped: list = []
        for arg in args:
            if isinstance(arg, Vector):
                unwrapped.extend(arg.components)
            else:
                unwrapped.append(arg)
        args = tuple(unwrapped)

        if all(arg is not None for arg in args):
            return self.function(*args)
        if all(arg is None for arg in args):
            raise ValueError("At least one argument must be non-None.")

        # Build a partially-applied function
        original = list(args)

        def partial(*z):
            filled = list(z)
            arguments = [filled.pop(0) if arg is None else arg for arg in original]
            return self(*arguments)

        n_free = sum(arg is None for arg in args)
        if n_free == 1:
            return Function(partial)
        return MultiVariableFunction(partial)


class BivariateFunction(MultiVariableFunction):
    """A :class:`MultiVariableFunction` specialised for exactly two arguments.

    This is a convenience alias — no additional behaviour beyond
    :class:`MultiVariableFunction`.
    """
