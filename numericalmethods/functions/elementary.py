from __future__ import annotations

import math
from typing import TYPE_CHECKING

from ..enums import InterpolationMethod, InterpolationForm
from .base import Function

if TYPE_CHECKING:
    pass


class Polynomial(Function):
    """A polynomial function :math:`p(x) = a_0 + a_1 x + \\cdots + a_n x^n`.

    Parameters
    ----------
    *coefficients:
        Coefficients in **ascending** degree order: ``a_0, a_1, ..., a_n``.
    """

    def __init__(self, *coefficients: float) -> None:
        self.coefficients = coefficients
        self.function = lambda x: sum(a * x**i for i, a in enumerate(coefficients))

    @staticmethod
    def interpolate(
        data: list[tuple[float, float]] | list[float],
        method: InterpolationMethod = InterpolationMethod.NEWTON,
        f: Function | None = None,
        form: InterpolationForm = InterpolationForm.STANDARD,
    ) -> Polynomial:
        """Construct an interpolating polynomial.

        Parameters
        ----------
        data:
            Either a list of ``(x, y)`` tuples, **or** a list of *x* values
            when *f* is provided.
        method:
            Interpolation algorithm (:class:`InterpolationMethod`).
        f:
            If provided, evaluates *f* at each *x* in *data* to build the
            ``(x, y)`` pairs.
        form:
            Newton form variant (:class:`InterpolationForm`): ``STANDARD``
            (divided differences), ``FORWARD_DIFF``, or ``BACKWARD_DIFF``.
        """
        if f is not None:
            data = [(x, f(x)) for x in data]

        if method is InterpolationMethod.LAGRANGE:
            return Polynomial._interpolate_lagrange(data)
        if method is InterpolationMethod.NEWTON:
            if form is InterpolationForm.STANDARD:
                return Polynomial._interpolate_newton(data)
            if form is InterpolationForm.FORWARD_DIFF:
                return Polynomial._interpolate_newton_forward_diff(data)
            if form is InterpolationForm.BACKWARD_DIFF:
                return Polynomial._interpolate_newton_backward_diff(data)

        raise ValueError(f"Unknown interpolation method/form combination: {method!r}/{form!r}")

    @staticmethod
    def _interpolate_lagrange(data: list[tuple[float, float]]) -> Polynomial:
        n = len(data)
        x = [d[0] for d in data]
        y = [d[1] for d in data]
        p = Polynomial(0)

        for i in range(n):
            L = Polynomial(1)
            for j in range(n):
                if i != j:
                    L = L * Polynomial(-x[j], 1) / Polynomial(x[i] - x[j])
            p = p + y[i] * L

        return p

    @staticmethod
    def _interpolate_newton(data: list[tuple[float, float]]) -> Polynomial:
        n = len(data)
        x = [d[0] for d in data]
        y = [d[1] for d in data]

        def divided_difference(i: int, j: int) -> float:
            if i == j:
                return y[i]
            return (divided_difference(i + 1, j) - divided_difference(i, j - 1)) / (x[j] - x[i])

        def factor_product(roots: list[float]) -> Polynomial:
            if not roots:
                return Polynomial(1)
            return Polynomial(-roots[0], 1) * factor_product(roots[1:])

        coefficients = [divided_difference(0, i) for i in range(n)]
        p = Polynomial(coefficients[0])
        for i in range(1, n):
            p = p + coefficients[i] * factor_product(x[:i])

        return p

    @staticmethod
    def _interpolate_newton_forward_diff(data: list[tuple[float, float]]) -> Polynomial:
        from ..util import Util

        data = sorted(data, key=lambda d: d[0])
        diffs = sorted([data[i + 1][0] - data[i][0] for i in range(len(data) - 1)])
        for j in range(len(diffs) - 1):
            if abs(diffs[j + 1] - diffs[j]) >= 1e-6:
                raise ValueError("x values must be equally spaced for Newton forward-difference.")

        h = abs(diffs[0])
        n = len(data) - 1

        p = Polynomial(data[0][1])
        for k in range(1, n + 1):
            p = p + Util.delta(k, 0, data) * Util.choose(Polynomial(-data[0][0] / h, 1 / h), k)

        return p

    @staticmethod
    def _interpolate_newton_backward_diff(data: list[tuple[float, float]]) -> Polynomial:
        from ..util import Util

        data = sorted(data, key=lambda d: d[0])
        diffs = sorted([data[i + 1][0] - data[i][0] for i in range(len(data) - 1)])
        for j in range(len(diffs) - 1):
            if abs(diffs[j + 1] - diffs[j]) >= 1e-6:
                raise ValueError("x values must be equally spaced for Newton backward-difference.")

        h = abs(diffs[0])
        n = len(data) - 1

        p = Polynomial(data[n][1])
        for k in range(1, n + 1):
            p = p + ((-1) ** k * Util.downdelta(k, n, data)) * Util.choose(
                Polynomial(data[n][0] / h, -1 / h), k
            )

        return p


class Exponent(Function):
    """The exponential function :math:`b^{f(x)}`.

    Parameters
    ----------
    f:
        Exponent function.
    base:
        Base of the exponential (defaults to :math:`e`).
    """

    def __init__(self, f: Function, base: float = math.e) -> None:
        self.function = lambda x: base ** f(x)


class Sin(Function):
    """The sine function :math:`\\sin(f(x))`."""

    def __init__(self, f: Function) -> None:
        self.function = lambda x: math.sin(f(x))


class Cos(Function):
    """The cosine function :math:`\\cos(f(x))`."""

    def __init__(self, f: Function) -> None:
        self.function = lambda x: math.cos(f(x))


class Tan(Function):
    """The tangent function :math:`\\tan(f(x))`."""

    def __init__(self, f: Function) -> None:
        self.function = lambda x: math.tan(f(x))


class Log(Function):
    """The logarithm :math:`\\log_b(f(x))`.

    Parameters
    ----------
    f:
        Argument function.
    base:
        Logarithm base (defaults to :math:`e`, i.e. natural logarithm).
    """

    def __init__(self, f: Function, base: float = math.e) -> None:
        self.function = lambda x: math.log(f(x), base)
