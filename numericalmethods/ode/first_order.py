from __future__ import annotations

from fractions import Fraction

from ..enums import ODEMethod
from ..functions.base import Function
from ..functions.multivariate import BivariateFunction, MultiVariableFunction
from ..functions.elementary import Polynomial
from ..linalg.vector import Vector
from .base import LinearODE


class FirstOrderLinearODE(LinearODE):
    """Numerical solver for a first-order linear IVP :math:`y'(x) = f(x, y(x))`.

    Parameters
    ----------
    f:
        Right-hand side :math:`f(x, y)`, given as a :class:`BivariateFunction`
        (or a :class:`~numericalmethods.linalg.Vector` of
        :class:`~numericalmethods.functions.multivariate.MultiVariableFunction` objects
        for systems of equations).
    a, b:
        Integration interval.
    y0:
        Initial condition :math:`y(a)`.
    """

    def __init__(
        self,
        f: BivariateFunction | Vector,
        a: float,
        b: float,
        y0: float | Vector,
    ) -> None:
        self.f = f
        self.a = a
        self.b = b
        self.y0 = y0

    def solve(
        self,
        h: float = 0.1,
        method: ODEMethod = ODEMethod.EULER,
        n: int = 1,
        step: int = 2,
        points: list[float] | None = None,
    ) -> Polynomial:
        """Solve the IVP and return an interpolating :class:`Polynomial`.

        Parameters
        ----------
        h:
            Step size.
        method:
            Solver algorithm (:class:`ODEMethod`).
        n:
            Order parameter for Taylor / Runge-Kutta methods.
        step:
            Number of steps for Adams methods.
        points:
            Known function values at previous grid points for multistep methods.
        """
        if points is None:
            points = []


        if method is ODEMethod.EULER:
            return self._solve_taylor(h, 1)
        if method is ODEMethod.RUNGE_KUTTA:
            return self._solve_runge_kutta(h, n)
        if method is ODEMethod.TAYLOR:
            return self._solve_taylor(h, n)
        if method is ODEMethod.TRAPEZOIDAL:
            return self._solve_trapezoidal(h)
        if method is ODEMethod.ADAMS_BASHFORTH:
            return self._solve_adams_bashforth(h, step, points)
        if method is ODEMethod.ADAMS_MOULTON:
            return self._solve_adams_moulton(h, step, points)
        if method is ODEMethod.PREDICTOR_CORRECTOR:
            return self._solve_predictor_corrector(h)

        raise ValueError(f"Unknown ODE method: {method!r}")

    # ------------------------------------------------------------------
    # Private solvers
    # ------------------------------------------------------------------

    def _solve_runge_kutta(self, h: float, n: int) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)

        if n == 1:
            return self.solve(h, method=ODEMethod.EULER)
        elif n == 2:
            for i in range(N):
                xi = self.a + i * h
                w.append(
                    w[i]
                    + (h / 2) * self.f(xi, w[i])
                    + (h / 2) * self.f(xi + h, w[i] + h * self.f(xi, w[i]))
                )
        elif n == 3:
            for i in range(N):
                xi = self.a + i * h
                k1 = self.f(xi, w[i])
                k2 = self.f(xi + (h / 3), w[i] + (h / 3) * k1)
                k3 = self.f(xi + (2 / 3) * h, w[i] + (2 / 3) * h * k2)
                w.append(w[i] + (h / 4) * (k1 + 3 * k3))
        elif n == 4:
            for i in range(N):
                xi = self.a + i * h
                k1 = h * self.f(xi, w[i])
                k2 = h * self.f(xi + h / 2, w[i] + 0.5 * k1)
                k3 = h * self.f(xi + h / 2, w[i] + 0.5 * k2)
                k4 = h * self.f(xi + h, w[i] + k3)
                w.append(w[i] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4))
        else:
            raise NotImplementedError("Runge-Kutta is only implemented for orders 1–4.")

        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except Exception:
            return w  # fall back to list for vector-valued problems

    def _solve_taylor(self, h: float, n: int) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)

        if n == 1:
            for i in range(N):
                xi = self.a + i * h
                w.append(w[i] + h * self.f(xi, w[i]))
        elif n == 2:
            for i in range(N):
                xi = self.a + i * h
                g = (
                    self.f(None, w[i]).differentiate()(xi)
                    + self.f(xi, w[i]) * self.f(xi, None).differentiate()(w[i])
                )
                w.append(w[i] + h * self.f(xi, w[i]) + (h**2) * g / 2)
        else:
            raise NotImplementedError("Taylor method is only implemented for orders 1 and 2.")

        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except Exception:
            return w

    def _solve_trapezoidal(self, h: float) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)

        for i in range(N):
            xi = self.a + i * h
            g = Function(
                lambda x, _h=h, _xi=xi, _wi=w[i]: _wi
                + (_h / 2) * (self.f(_xi, _wi) + self.f(_xi + _h, x))
            )
            w.append(g.fixed_point(w[i]))

        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except Exception:
            return w

    @staticmethod
    def _adams_bashforth_coefficients(k: int) -> list[Fraction]:
        """Compute Adams-Bashforth *k*-step coefficients via Lagrange interpolation.

        Returns ``[b_0, ..., b_{k-1}]`` (as :class:`~fractions.Fraction` objects) such
        that :math:`w_{n+1} = w_n + h \\sum_{j=0}^{k-1} b_j f(t_{n-j}, w_{n-j})`.
        """
        nodes = [Fraction(-j) for j in range(k)]
        coefficients = []
        for m in range(k):
            poly = [Fraction(1)]
            for j in range(k):
                if j != m:
                    new_poly = [Fraction(0)] * (len(poly) + 1)
                    for idx, c in enumerate(poly):
                        new_poly[idx + 1] += c
                        new_poly[idx] -= nodes[j] * c
                    poly = new_poly
            denom = Fraction(1)
            for j in range(k):
                if j != m:
                    denom *= nodes[m] - nodes[j]
            integral = sum(c / (i + 1) for i, c in enumerate(poly))
            coefficients.append(integral / denom)
        return coefficients

    @staticmethod
    def _adams_moulton_coefficients(k: int) -> list[Fraction]:
        """Compute Adams-Moulton *k*-step coefficients via Lagrange interpolation.

        Returns ``[b_0, ..., b_k]`` (as :class:`~fractions.Fraction` objects) such that
        :math:`w_{n+1} = w_n + h(b_0 f(t_{n+1}, w_{n+1}) + \\sum_{j=1}^{k} b_j f(t_{n+1-j}, w_{n+1-j}))`.
        """
        nodes = [Fraction(1 - j) for j in range(k + 1)]
        coefficients = []
        for m in range(k + 1):
            poly = [Fraction(1)]
            for j in range(k + 1):
                if j != m:
                    new_poly = [Fraction(0)] * (len(poly) + 1)
                    for idx, c in enumerate(poly):
                        new_poly[idx + 1] += c
                        new_poly[idx] -= nodes[j] * c
                    poly = new_poly
            denom = Fraction(1)
            for j in range(k + 1):
                if j != m:
                    denom *= nodes[m] - nodes[j]
            integral = sum(c / (i + 1) for i, c in enumerate(poly))
            coefficients.append(integral / denom)
        return coefficients

    def _solve_adams_bashforth(
        self, h: float, step: int, points: list[float]
    ) -> Polynomial:
        if step < 2:
            raise ValueError("Adams-Bashforth requires at least 2 steps.")
        w = [self.y0] + list(points)
        N = int((self.b - self.a) / h)
        coeffs = [float(c) for c in FirstOrderLinearODE._adams_bashforth_coefficients(step)]

        for i in range(step - 1, N):
            xi = self.a + i * h
            w.append(
                w[i]
                + h * sum(coeffs[j] * self.f(xi - j * h, w[i - j]) for j in range(step))
            )

        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except Exception:
            return w

    def _solve_adams_moulton(
        self, h: float, step: int, points: list[float]
    ) -> Polynomial:
        if step < 2:
            raise ValueError("Adams-Moulton requires at least 2 steps.")
        w = [self.y0] + list(points)
        N = int((self.b - self.a) / h)
        coeffs = [float(c) for c in FirstOrderLinearODE._adams_moulton_coefficients(step)]

        for i in range(step - 1, N):
            xi = self.a + i * h
            explicit = sum(
                coeffs[j] * self.f(xi - (j - 1) * h, w[i - j + 1]) for j in range(1, step + 1)
            )
            c0 = coeffs[0]
            g = Function(
                lambda x, _h=h, _xi=xi, _c0=c0, _explicit=explicit, _wi=w[i]: (
                    _wi + _h * (_c0 * self.f(_xi + _h, x) + _explicit)
                )
            )
            w.append(g.fixed_point(w[i]))

        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except Exception:
            return w

    def _solve_predictor_corrector(self, h: float) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)
        alphas = [self.a + h * i for i in range(1, 4)]

        # Starting values via RK4
        ivp = FirstOrderLinearODE(self.f, self.a, self.a + 3 * h, self.y0)
        rk_sol = ivp._solve_runge_kutta(h, 4)
        for xi in alphas:
            w.append(rk_sol(xi))

        for i in range(3, N):
            xi = self.a + i * h
            prediction = w[i] + (h / 24) * (
                55 * self.f(xi, w[i])
                - 59 * self.f(xi - h, w[i - 1])
                + 37 * self.f(xi - 2 * h, w[i - 2])
                - 9 * self.f(xi - 3 * h, w[i - 3])
            )
            correction = w[i] + (h / 24) * (
                9 * self.f(xi + h, prediction)
                + 19 * self.f(xi, w[i])
                - 5 * self.f(xi - h, w[i - 1])
                + self.f(xi - 2 * h, w[i - 2])
            )
            w.append(correction)

        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except Exception:
            return w
