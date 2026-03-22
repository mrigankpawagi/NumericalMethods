from __future__ import annotations

from ..enums import BVPMethod, NonlinearBVPMethod, ODEMethod
from ..functions.base import Function
from ..functions.elementary import Polynomial
from ..functions.multivariate import MultiVariableFunction
from ..linalg.vector import Vector
from ..linalg.matrix import Matrix
from ..linalg.linear_system import LinearSystem
from .base import LinearODE, OrdinaryDifferentialEquation
from .first_order import FirstOrderLinearODE


class SecondOrderLinearODE_BVP(LinearODE):
    """Solver for the second-order linear BVP

    .. math::

        y''(x) = p(x)\\,y'(x) + q(x)\\,y(x) + r(x),
        \\quad y(a) = y_0,\\quad y(b) = y_1.

    Parameters
    ----------
    p, q, r:
        Coefficient functions.
    a, b:
        Boundary points.
    y0, y1:
        Boundary values :math:`y(a)` and :math:`y(b)`.
    """

    def __init__(
        self,
        p: Function,
        q: Function,
        r: Function,
        a: float,
        b: float,
        y0: float,
        y1: float,
    ) -> None:
        self.p = p
        self.q = q
        self.r = r
        self.a = a
        self.b = b
        self.y0 = y0
        self.y1 = y1

    def solve(
        self,
        h: float = 0.1,
        method: BVPMethod = BVPMethod.SHOOTING,
    ) -> Polynomial:
        """Solve the BVP and return an interpolating :class:`Polynomial`.

        Parameters
        ----------
        h:
            Step size / mesh width.
        method:
            Algorithm (:class:`BVPMethod`).
        """

        if method is BVPMethod.SHOOTING:
            return self._solve_shooting(h)
        if method is BVPMethod.FINITE_DIFFERENCE:
            return self._solve_finite_difference(h)

        raise ValueError(f"Unknown BVP method: {method!r}")

    def _solve_shooting(self, h: float) -> Polynomial:
        IVP1 = SecondOrderODE_IVP(
            MultiVariableFunction(
                lambda t, u1, u2: self.p(t) * u2 + self.q(t) * u1 + self.r(t)
            ),
            self.a, self.b, self.y0, 0,
        )
        IVP2 = SecondOrderODE_IVP(
            MultiVariableFunction(lambda t, u1, u2: self.p(t) * u2 + self.q(t) * u1),
            self.a, self.b, 0, 1,
        )

        sol1 = IVP1.solve(h)
        sol2 = IVP2.solve(h)
        c = (self.y1 - sol1(self.b)) / sol2(self.b)

        return sol1 + c * sol2

    def _solve_finite_difference(self, h: float) -> Polynomial:
        N = int((self.b - self.a) / h) - 1
        A = Matrix(*[Vector(*[0.0] * N) for _ in range(N)])
        b = Vector(*[0.0] * N)

        # First row
        A[0][0] = -(2 + h**2 * self.q(self.a + h))
        A[0][1] = 1 - (h / 2) * self.p(self.a + h)
        b[0] = h**2 * self.r(self.a + h) - (1 + (h / 2) * self.p(self.a + h)) * self.y0

        # Middle rows
        for i in range(1, N - 1):
            xi = self.a + (i + 1) * h
            A[i][i - 1] = 1 + (h / 2) * self.p(xi)
            A[i][i] = -(2 + h**2 * self.q(xi))
            A[i][i + 1] = 1 - (h / 2) * self.p(xi)
            b[i] = h**2 * self.r(xi)

        # Last row
        A[N - 1][N - 2] = 1 + (h / 2) * self.p(self.b - h)
        A[N - 1][N - 1] = -(2 + h**2 * self.q(self.b - h))
        b[N - 1] = (
            h**2 * self.r(self.b - h)
            - (1 - (h / 2) * self.p(self.b - h)) * self.y1
        )

        sol = LinearSystem(A, b).solve()
        w = [self.y0] + sol.components + [self.y1]

        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 2)])


class SecondOrderODE_IVP(OrdinaryDifferentialEquation):
    """Solver for the second-order nonlinear IVP

    .. math::

        y''(x) = f(x,\\,y(x),\\,y'(x)),
        \\quad y(a) = y_0,\\quad y'(a) = y_1.

    The equation is reduced to a first-order system before numerical integration.

    Parameters
    ----------
    f:
        Right-hand side :math:`f(x, y, y')`.
    a, b:
        Integration interval.
    y0:
        Initial value :math:`y(a)`.
    y1:
        Initial slope :math:`y'(a)`.
    """

    def __init__(
        self,
        f: MultiVariableFunction,
        a: float,
        b: float,
        y0: float,
        y1: float,
    ) -> None:
        self.f = f
        self.a = a
        self.b = b
        self.y0 = y0
        self.y1 = y1

    def solve(
        self,
        h: float = 0.1,
        method: ODEMethod = ODEMethod.EULER,
        n: int = 1,
        step: int = 2,
        points: list[float] | None = None,
    ) -> Polynomial:
        """Solve the IVP via system reduction and return an interpolating :class:`Polynomial`.

        Accepts the same *method*, *n*, *step*, and *points* parameters as
        :meth:`FirstOrderLinearODE.solve`.
        """
        if points is None:
            points = []

        U0 = Vector(self.y0, self.y1)
        F = Vector(
            MultiVariableFunction(lambda t, u1, u2: u2),
            MultiVariableFunction(lambda t, u1, u2: self.f(t, u1, u2)),
        )
        ivp = FirstOrderLinearODE(F, self.a, self.b, U0)
        sol = ivp.solve(h, method, n, step, points)

        w = [x[0] for x in sol]
        N = int((self.b - self.a) / h)
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(len(w))])


class SecondOrderODE_BVP(OrdinaryDifferentialEquation):
    """Solver for the second-order nonlinear BVP

    .. math::

        y''(x) = f(x,\\,y(x),\\,y'(x)),
        \\quad y(a) = y_0,\\quad y(b) = y_1.

    Parameters
    ----------
    f:
        Right-hand side :math:`f(x, y, y')`.
    a, b:
        Boundary points.
    y0, y1:
        Boundary values :math:`y(a)` and :math:`y(b)`.
    """

    def __init__(
        self,
        f: MultiVariableFunction,
        a: float,
        b: float,
        y0: float,
        y1: float,
    ) -> None:
        self.f = f
        self.a = a
        self.b = b
        self.y0 = y0
        self.y1 = y1

    def solve(
        self,
        h: float = 0.1,
        method: NonlinearBVPMethod = NonlinearBVPMethod.SHOOTING_NEWTON,
        M: int = 100,
        tol: float = 1e-5,
        initial_approximation=None,
    ) -> Polynomial | None:
        """Solve the nonlinear BVP and return an interpolating :class:`Polynomial`.

        Parameters
        ----------
        h:
            Step size / mesh width.
        method:
            Algorithm (:class:`NonlinearBVPMethod`).
        M:
            Maximum number of outer iterations.
        tol:
            Convergence tolerance.
        initial_approximation:
            Initial guess for :math:`y'(a)` (shooting) or interior node values
            (finite difference).
        """

        if method is NonlinearBVPMethod.SHOOTING_NEWTON:
            return self._solve_shooting_newton(h, M, tol, initial_approximation)
        if method is NonlinearBVPMethod.FINITE_DIFFERENCE:
            return self._solve_finite_difference(h, M, tol)

        raise ValueError(f"Unknown nonlinear BVP method: {method!r}")

    def _solve_shooting_newton(
        self, h: float, M: int, tol: float, initial_approximation
    ) -> Polynomial | None:
        t = 1.0 if initial_approximation is None else initial_approximation
        for _ in range(M):
            IVP1 = SecondOrderODE_IVP(
                MultiVariableFunction(lambda t_var, u1, u2: self.f(t_var, u1, u2)),
                self.a, self.b, self.y0, t,
            )
            y = IVP1.solve(h)

            p = Function(
                lambda x: self.f(x, None, y.differentiate()(x)).differentiate()(y(x))
            )
            q = Function(
                lambda x: self.f(x, y(x), None).differentiate()(y.differentiate()(x))
            )
            r = Function(lambda x: 0)

            IVP2 = SecondOrderLinearODE_BVP(p, q, r, self.a, self.b, 0, 1)
            z = IVP2.solve(h)

            t_new = t - (y(self.b) - self.y1) / z(self.b)
            if abs(t_new - t) < tol:
                return y
            t = t_new

        return None

    def _solve_finite_difference(
        self, h: float, M: int, tol: float
    ) -> Polynomial | None:
        N = int(round((self.b - self.a) / h)) - 1
        pd_step = 1e-5

        def df_dy(x: float, y: float, z: float) -> float:
            return (self.f(x, y + pd_step, z) - self.f(x, y - pd_step, z)) / (2 * pd_step)

        def df_dz(x: float, y: float, z: float) -> float:
            return (self.f(x, y, z + pd_step) - self.f(x, y, z - pd_step)) / (2 * pd_step)

        # Linear interpolation as initial guess
        u = [
            self.y0 + (j + 1) * h * (self.y1 - self.y0) / (self.b - self.a)
            for j in range(N)
        ]

        for _ in range(M):
            A = Matrix(*[Vector(*[0.0] * N) for _ in range(N)])
            F = [0.0] * N

            for j in range(N):
                xi = self.a + (j + 1) * h
                wi = u[j]
                wi_prev = u[j - 1] if j > 0 else self.y0
                wi_next = u[j + 1] if j < N - 1 else self.y1
                zprime = (wi_next - wi_prev) / (2 * h)

                F[j] = wi_next - 2 * wi + wi_prev - h**2 * self.f(xi, wi, zprime)
                A[j][j] = -2 - h**2 * df_dy(xi, wi, zprime)
                if j > 0:
                    A[j][j - 1] = 1 + (h / 2) * df_dz(xi, wi, zprime)
                if j < N - 1:
                    A[j][j + 1] = 1 - (h / 2) * df_dz(xi, wi, zprime)

            delta = LinearSystem(A, Vector(*[-fi for fi in F])).solve()
            u = [u[j] + delta[j] for j in range(N)]

            if max(abs(delta[j]) for j in range(N)) < tol:
                w = [self.y0] + u + [self.y1]
                return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 2)])

        return None
