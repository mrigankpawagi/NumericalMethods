from __future__ import annotations

from ..enums import LinearSolverMethod
from .vector import Vector
from .matrix import Matrix


class LinearSystem:
    """A system of linear equations :math:`Ax = b`.

    Parameters
    ----------
    A:
        The coefficient :class:`Matrix`.
    b:
        The right-hand-side :class:`Vector`.
    """

    def __init__(self, A: Matrix, b: Vector) -> None:
        self.A = A
        self.b = b
        self.N = len(A)

        if len(b) != self.N:
            raise ValueError("A and b must have the same number of rows.")

    def solve(
        self,
        method: LinearSolverMethod | str = LinearSolverMethod.GAUSS_ELIMINATION,
        tol: float = 1e-5,
        initial_approximation: Vector | None = None,
        max_iterations: int = 100,
    ) -> Vector:
        """Solve the linear system and return the solution :class:`Vector`.

        Parameters
        ----------
        method:
            Solver algorithm to use.  Accepts :class:`LinearSolverMethod` enum
            values or their string equivalents (``"gauss_elimination"``,
            ``"gauss_jacobi"``, ``"gauss_seidel"``).
        tol:
            Convergence tolerance for iterative solvers.
        initial_approximation:
            Starting guess for iterative solvers.
        max_iterations:
            Maximum number of iterations for iterative solvers.
        """
        method = LinearSolverMethod(method)

        if method is LinearSolverMethod.GAUSS_ELIMINATION:
            return self._gauss_elimination()
        if method is LinearSolverMethod.GAUSS_JACOBI:
            return self._gauss_jacobi(tol, initial_approximation, max_iterations)
        if method is LinearSolverMethod.GAUSS_SEIDEL:
            return self._gauss_seidel(tol, initial_approximation, max_iterations)

        raise ValueError(f"Unknown method: {method!r}")

    # ------------------------------------------------------------------
    # Private solvers
    # ------------------------------------------------------------------

    def _gauss_elimination(self) -> Vector:
        for i in range(self.N - 1):
            p = next((j for j in range(i, self.N) if abs(self.A[j][i]) != 0), None)
            if p is None:
                raise ValueError("No unique solution exists.")

            if p != i:
                self.A[i], self.A[p] = self.A[p], self.A[i]
                self.b[i], self.b[p] = self.b[p], self.b[i]

            for j in range(i + 1, self.N):
                m = self.A[j][i] / self.A[i][i]
                self.A[j] = self.A[j] - m * self.A[i]
                self.b[j] = self.b[j] - m * self.b[i]

        if abs(self.A[self.N - 1][self.N - 1]) == 0:
            raise ValueError("No unique solution exists.")

        x: list[float] = [0.0] * self.N
        x[self.N - 1] = self.b[self.N - 1] / self.A[self.N - 1][self.N - 1]
        for i in range(self.N - 2, -1, -1):
            x[i] = (
                self.b[i] - sum(self.A[i][j] * x[j] for j in range(i + 1, self.N))
            ) / self.A[i][i]

        return Vector(*x)

    def _gauss_jacobi(
        self, tol: float, initial_approximation: Vector | None, max_iterations: int
    ) -> Vector | None:
        if initial_approximation is None:
            raise ValueError("initial_approximation must be provided for iterative solvers.")

        x0 = initial_approximation
        for _ in range(max_iterations):
            x1 = [
                (self.b[i] - sum(self.A[i][j] * x0[j] for j in range(self.N) if j != i))
                / self.A[i][i]
                for i in range(self.N)
            ]
            if (
                max(abs(x1[i] - x0[i]) for i in range(self.N))
                / max(abs(x1[i]) for i in range(self.N))
                < tol
            ):
                return Vector(*x1)
            x0 = x1

        return None

    def _gauss_seidel(
        self, tol: float, initial_approximation: Vector | None, max_iterations: int
    ) -> Vector | None:
        if initial_approximation is None:
            raise ValueError("initial_approximation must be provided for iterative solvers.")

        x0 = initial_approximation
        for _ in range(max_iterations):
            x1 = [0.0] * self.N
            for i in range(self.N):
                x1[i] = (
                    self.b[i]
                    - sum(self.A[i][j] * x1[j] for j in range(i))
                    - sum(self.A[i][j] * x0[j] for j in range(i + 1, self.N))
                ) / self.A[i][i]

            if (
                max(abs(x1[i] - x0[i]) for i in range(self.N))
                / max(abs(x1[i]) for i in range(self.N))
                < tol
            ):
                return Vector(*x1)
            x0 = x1

        return None
