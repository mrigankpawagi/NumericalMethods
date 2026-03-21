from __future__ import annotations

import math
from fractions import Fraction
from typing import TYPE_CHECKING, Callable

from ..enums import DifferentiationMethod, IntegrationMethod, RootFindingMethod

if TYPE_CHECKING:
    pass


class Function:
    """A callable mathematical function with composable arithmetic and numerical methods.

    A :class:`Function` wraps any callable ``f: float -> float`` and provides:

    * **Arithmetic operators** — ``+``, ``-``, ``*``, ``/``, ``**`` (and their
      reflected variants) that return new :class:`Function` objects.
    * **Differentiation** — :meth:`differentiate` and :meth:`multi_differentiate`
      via finite differences, or by supplying an exact derivative.
    * **Integration** — :meth:`integrate` with several quadrature rules.
    * **Root finding** — :meth:`root` with bisection, Newton-Raphson, secant,
      Regula Falsi, and modified Newton methods.
    * **Fixed-point iteration** — :meth:`fixed_point`.
    * **Plotting** — :meth:`plot` (requires ``matplotlib``).

    Parameters
    ----------
    function:
        Any callable ``f(x) -> float``.
    """

    def __init__(self, function: Callable) -> None:
        self.function = function

    # ------------------------------------------------------------------
    # Calling
    # ------------------------------------------------------------------

    def __call__(self, x):
        if callable(x):
            return Function(lambda p: self(x(p)))
        return self.function(x)

    # ------------------------------------------------------------------
    # Arithmetic operators
    # ------------------------------------------------------------------

    def __add__(self, other: Function) -> Function:
        return Function(lambda x: self(x) + other(x))

    def __sub__(self, other: Function) -> Function:
        return Function(lambda x: self(x) - other(x))

    def __mul__(self, other: Function) -> Function:
        return Function(lambda x: self(x) * other(x))

    def __truediv__(self, other: Function) -> Function:
        return Function(lambda x: self(x) / other(x))

    def __pow__(self, other: Function) -> Function:
        return Function(lambda x: self(x) ** other(x))

    def __radd__(self, other: float) -> Function:
        return Function(lambda x: other + self(x))

    def __rsub__(self, other: float) -> Function:
        return Function(lambda x: other - self(x))

    def __rmul__(self, other: float) -> Function:
        return Function(lambda x: other * self(x))

    def __rtruediv__(self, other: float) -> Function:
        return Function(lambda x: other / self(x))

    def __neg__(self) -> Function:
        return Function(lambda x: -self(x))

    # ------------------------------------------------------------------
    # Differentiation
    # ------------------------------------------------------------------

    def differentiate(
        self,
        func: Function | Callable | None = None,
        h: float = 1e-5,
        method: DifferentiationMethod | str = DifferentiationMethod.FORWARD,
    ) -> Function | None:
        """Return a numerical derivative, or set an exact derivative.

        * If *func* is ``None`` (default), returns a :class:`Function`
          approximating the derivative via finite differences.
        * If *func* is a :class:`Function` or callable, stores it as the
          exact derivative (accessible in Newton-type root finders).

        Parameters
        ----------
        func:
            Exact derivative function, or ``None`` to use finite differences.
        h:
            Step size for finite-difference approximation.
        method:
            Finite-difference scheme: ``FORWARD``, ``BACKWARD``, or ``CENTRAL``.
        """
        if isinstance(func, Function):
            self.derivative = func
            return None
        if callable(func):
            self.derivative = Function(func)
            return None

        method = DifferentiationMethod(method)
        if method is DifferentiationMethod.FORWARD:
            return self._differentiate_forward(h)
        if method is DifferentiationMethod.BACKWARD:
            return self._differentiate_forward(-h)
        if method is DifferentiationMethod.CENTRAL:
            return self._differentiate_central(h)

        raise ValueError(f"Unknown differentiation method: {method!r}")

    def _differentiate_forward(self, h: float) -> Function:
        return Function(lambda x: (self(x + h) - self(x)) / h)

    def _differentiate_central(self, h: float) -> Function:
        return Function(lambda x: (self(x + h) - self(x - h)) / (2 * h))

    def multi_differentiate(
        self,
        n: int,
        h: float = 1e-5,
        method: DifferentiationMethod | str = DifferentiationMethod.FORWARD,
    ) -> Function:
        """Return the *n*-th order derivative via repeated finite differences.

        Parameters
        ----------
        n:
            Order of differentiation.
        h:
            Step size.
        method:
            Finite-difference scheme.
        """
        if n == 0:
            return self
        return self.differentiate(h=h, method=method).multi_differentiate(n - 1, h, method)

    # ------------------------------------------------------------------
    # Integration
    # ------------------------------------------------------------------

    def integral(self, func: Function | Callable | None = None) -> None:
        """Set an exact antiderivative.

        Parameters
        ----------
        func:
            A :class:`Function` (or callable) representing the antiderivative.
            Once set, :meth:`integrate` with no *method* will use it.
        """
        if isinstance(func, Function):
            self._integral = func
        elif callable(func):
            self._integral = Function(func)
        else:
            raise NotImplementedError("Automatic antiderivative computation is not implemented.")

    def integrate(
        self,
        a: float,
        b: float,
        method: IntegrationMethod | str | None = None,
        n: int | None = None,
    ) -> float:
        """Compute the definite integral :math:`\\int_a^b f(x)\\,dx`.

        Parameters
        ----------
        a, b:
            Integration bounds.
        method:
            Quadrature rule.  Accepts :class:`IntegrationMethod` values or
            strings: ``"rectangular"``, ``"midpoint"``, ``"trapezoidal"``,
            ``"simpson"``, ``"gauss"``.  If ``None`` and an exact antiderivative
            has been set via :meth:`integral`, it is used instead.
        n:
            Number of subintervals (or quadrature points for Gauss).
        """
        if method is not None:
            method = IntegrationMethod(method)

        if method is IntegrationMethod.RECTANGULAR:
            return self._integrate_rectangular(a, b, n)
        if method is IntegrationMethod.MIDPOINT:
            return self._integrate_midpoint(a, b, n)
        if method is IntegrationMethod.TRAPEZOIDAL:
            return self._integrate_trapezoidal(a, b, n)
        if method is IntegrationMethod.SIMPSON:
            return self._integrate_simpson(a, b, n)
        if method is IntegrationMethod.GAUSS:
            return self._integrate_gauss(a, b, n)

        if hasattr(self, "_integral"):
            return self._integral(b) - self._integral(a)

        raise ValueError("Specify a method or set an exact antiderivative via .integral().")

    def _integrate_rectangular(self, a: float, b: float, n: int | None) -> float:
        if not n:
            return (b - a) * self(a)
        h = (b - a) / n
        return h * sum(self(a + i * h) for i in range(n))

    def _integrate_midpoint(self, a: float, b: float, n: int | None) -> float:
        if not n:
            return (b - a) * self((a + b) / 2)
        h = (b - a) / n
        return h * sum(self(a + i * h + h / 2) for i in range(n))

    def _integrate_trapezoidal(self, a: float, b: float, n: int | None) -> float:
        if not n:
            return (b - a) * (self(a) + self(b)) / 2
        h = (b - a) / n
        return h * (self(a) + 2 * sum(self(a + i * h) for i in range(1, n)) + self(b)) / 2

    def _integrate_simpson(self, a: float, b: float, n: int | None) -> float:
        if not n:
            return (b - a) * (self(a) + 4 * self((a + b) / 2) + self(b)) / 6
        h = (b - a) / n
        return (
            h
            * (
                self(a)
                + 4 * sum(self(a + i * h + h / 2) for i in range(n))
                + 2 * sum(self(a + i * h) for i in range(1, n))
                + self(b)
            )
            / 6
        )

    def _integrate_gauss(self, a: float, b: float, n: int | None) -> float:
        # Import here to avoid circular imports at module level
        from .elementary import Polynomial

        t = Polynomial((a + b) / 2, (b - a) / 2)
        g = ((b - a) / 2) * self(t)
        if n == 1:
            return 2 * g(0)
        if n == 2:
            return g(-1 / math.sqrt(3)) + g(1 / math.sqrt(3))

        raise NotImplementedError("Gaussian quadrature is only implemented for n=1 and n=2.")

    # ------------------------------------------------------------------
    # Root finding
    # ------------------------------------------------------------------

    def root(
        self,
        method: RootFindingMethod | str,
        a: float | None = None,
        b: float | None = None,
        p0: float | None = None,
        p1: float | None = None,
        TOLERANCE: float = 1e-10,
        N: int = 100,
        return_iterations: bool = False,
        early_stop: int | None = None,
    ) -> float | tuple[float, int] | None:
        """Find a root of the function.

        Parameters
        ----------
        method:
            Algorithm to use.  Accepts :class:`RootFindingMethod` values or
            strings: ``"bisection"``, ``"newton"``, ``"modified_newton"``,
            ``"secant"``, ``"regula_falsi"``.
        a, b:
            Bracket for bisection (both required).
        p0, p1:
            Initial guesses (requirements vary by method).
        TOLERANCE:
            Convergence tolerance.
        N:
            Maximum number of iterations.
        return_iterations:
            If ``True``, returns ``(root, n_iterations)`` instead of just the root.
        early_stop:
            Stop after this many iterations regardless of convergence (useful for
            computing the *k*-th iterate explicitly).
        """
        method = RootFindingMethod(method)

        if method is RootFindingMethod.BISECTION:
            if a is None or b is None:
                raise ValueError("Bisection requires both a and b.")
            if a >= b:
                raise ValueError("a must be less than b.")
            if self(a) * self(b) >= 0:
                raise ValueError("f(a) and f(b) must have opposite signs.")
            sol, n = self._bisection(a, b, TOLERANCE, N, early_stop)
        elif method is RootFindingMethod.NEWTON:
            if p0 is None:
                raise ValueError("Newton's method requires p0.")
            sol, n = self._newton(p0, TOLERANCE, N, early_stop)
        elif method is RootFindingMethod.SECANT:
            if p0 is None or p1 is None:
                raise ValueError("Secant method requires both p0 and p1.")
            sol, n = self._secant(p0, p1, TOLERANCE, N, early_stop)
        elif method is RootFindingMethod.REGULA_FALSI:
            if p0 is None or p1 is None:
                raise ValueError("Regula falsi requires both p0 and p1.")
            if self(p0) * self(p1) >= 0:
                raise ValueError("f(p0) and f(p1) must have opposite signs.")
            sol, n = self._regula_falsi(p0, p1, TOLERANCE, N, early_stop)
        elif method is RootFindingMethod.MODIFIED_NEWTON:
            if p0 is None:
                raise ValueError("Modified Newton's method requires p0.")
            sol, n = self._modified_newton(p0, TOLERANCE, N, early_stop)
        else:
            raise ValueError(f"Unknown root-finding method: {method!r}")

        if return_iterations:
            return sol, n
        return sol

    def _bisection(
        self, a: float, b: float, tol: float, N: int, early_stop: int | None
    ) -> tuple[float | None, int]:
        for i in range(N):
            p = (a + b) / 2
            if self(p) == 0 or abs(a - b) < tol or (early_stop is not None and i >= early_stop):
                return p, i + 1
            if self(a) * self(p) > 0:
                a = p
            else:
                b = p
        return None, N

    def _newton(
        self, p0: float, tol: float, N: int, early_stop: int | None
    ) -> tuple[float | None, int]:
        deriv = self.differentiate()
        try:
            for i in range(N):
                p = p0 - self(p0) / deriv(p0)
                if abs(p - p0) < tol or (early_stop is not None and i >= early_stop):
                    return p, i + 1
                p0 = p
        except (ZeroDivisionError, OverflowError):
            return None, i  # noqa: F821  (i is defined when exception fires)
        return None, N

    def _modified_newton(
        self, p0: float, tol: float, N: int, early_stop: int | None
    ) -> tuple[float | None, int]:
        deriv = self.differentiate()
        double_deriv = deriv.differentiate()
        try:
            for i in range(N):
                p = p0 - self(p0) * deriv(p0) / (deriv(p0) ** 2 - self(p0) * double_deriv(p0))
                if abs(p - p0) < tol or (early_stop is not None and i >= early_stop):
                    return p, i + 1
                p0 = p
        except (ZeroDivisionError, OverflowError):
            return None, i  # noqa: F821
        return None, N

    def _secant(
        self, p0: float, p1: float, tol: float, N: int, early_stop: int | None
    ) -> tuple[float | None, int]:
        for i in range(N):
            p = p1 - self(p1) * (p1 - p0) / (self(p1) - self(p0))
            if abs(p - p1) < tol or (early_stop is not None and i >= early_stop):
                return p, i + 1
            p0 = p1
            p1 = p
        return None, N

    def _regula_falsi(
        self, p0: float, p1: float, tol: float, N: int, early_stop: int | None
    ) -> tuple[float | None, int]:
        for i in range(N):
            p = p1 - self(p1) * (p1 - p0) / (self(p1) - self(p0))
            if abs(p - p1) < tol or (early_stop is not None and i >= early_stop):
                return p, i + 1
            if self(p0) * self(p) > 0:
                p0 = p1
            p1 = p
        return None, N

    def fixed_point(
        self, p0: float, TOLERANCE: float = 1e-10, N: int = 100
    ) -> float | None:
        """Apply fixed-point iteration :math:`p_{n+1} = f(p_n)`.

        Parameters
        ----------
        p0:
            Initial approximation.
        TOLERANCE:
            Convergence tolerance.
        N:
            Maximum number of iterations.
        """
        try:
            for _ in range(N):
                p = self(p0)
                if abs(p - p0) < TOLERANCE:
                    return p
                p0 = p
        except OverflowError:
            return None
        return None

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def plot(
        self, min: float, max: float, N: int = 1000, file: str = "", clear: bool = False
    ) -> None:
        """Plot the function on ``[min, max]``.

        Requires ``matplotlib``.

        Parameters
        ----------
        min, max:
            Domain bounds.
        N:
            Number of sample points.
        file:
            If non-empty, save to this path instead of displaying.
        clear:
            If ``True``, clear the current figure before plotting.
        """
        import matplotlib.pyplot as plt

        x = [min + (i / N) * (max - min) for i in range(N)]
        y = [self(t) for t in x]

        if clear:
            plt.clf()
        plt.plot(x, y)
        plt.xlabel("x")
        plt.ylabel("y")
        if file:
            plt.savefig(file)
        else:
            plt.show()
