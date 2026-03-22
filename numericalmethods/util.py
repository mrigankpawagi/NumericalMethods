from __future__ import annotations

import math

from .functions.base import Function
from .functions.elementary import Polynomial


class Util:
    """Utility helpers for interpolation arithmetic."""

    @staticmethod
    def choose(s: Function, k: int) -> Polynomial:
        """Return the falling-factorial polynomial :math:`\\binom{s}{k}`.

        Parameters
        ----------
        s:
            A :class:`~numericalmethods.functions.base.Function` (typically a
            linear polynomial in *x*).
        k:
            Non-negative integer.
        """
        res = Polynomial(1)
        for i in range(k):
            res = res * (-i + s)
        return (1 / math.factorial(k)) * res

    @staticmethod
    def delta(k: int, i: int, data: list[tuple[float, float]]) -> float:
        """Compute the *k*-th forward difference :math:`\\Delta^k y_i`.

        Parameters
        ----------
        k:
            Order of the forward difference.
        i:
            Starting index.
        data:
            Sequence of ``(x, y)`` pairs.
        """
        if k == 1:
            return data[i + 1][1] - data[i][1]
        return Util.delta(k - 1, i + 1, data) - Util.delta(k - 1, i, data)

    @staticmethod
    def downdelta(k: int, i: int, data: list[tuple[float, float]]) -> float:
        """Compute the *k*-th backward difference :math:`\\nabla^k y_i`.

        Parameters
        ----------
        k:
            Order of the backward difference.
        i:
            Starting index.
        data:
            Sequence of ``(x, y)`` pairs.
        """
        if k == 1:
            return data[i][1] - data[i - 1][1]
        return Util.downdelta(k - 1, i, data) - Util.downdelta(k - 1, i - 1, data)
    