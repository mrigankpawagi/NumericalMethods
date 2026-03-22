from __future__ import annotations

from ..linalg.linear_system import LinearSystem
from ..linalg.matrix import Matrix
from ..linalg.vector import Vector
from ..enums import LinearSolverMethod
from .multivariate import BivariateFunction


class FEM2D:
    r"""Solve 2D Poisson problems on rectangles with linear triangular FEM.

    Solves

    .. math::
        -\Delta u(x,y) = f(x,y)\ 	ext{in}\ \Omega=[x_a,x_b]	imes[y_a,y_b],
        \quad u=g\ 	ext{on}\ \partial\Omega.

    Parameters
    ----------
    f:
        Source term :math:`f(x,y)`.
    g:
        Dirichlet boundary condition :math:`g(x,y)`.
    xa, xb, ya, yb:
        Rectangle bounds.
    """

    def __init__(self, f: BivariateFunction, g: BivariateFunction, xa: float, xb: float, ya: float, yb: float) -> None:
        self.f = f
        self.g = g
        self.xa = xa
        self.xb = xb
        self.ya = ya
        self.yb = yb

    def solve(self, nx: int = 4, ny: int = 4) -> BivariateFunction:
        """Return a :class:`BivariateFunction` FEM approximation on an ``nx x ny`` mesh."""
        hx = (self.xb - self.xa) / nx
        hy = (self.yb - self.ya) / ny
        n_total = (nx + 1) * (ny + 1)

        def nidx(i: int, j: int) -> int:
            return i * (ny + 1) + j

        def ncoords(i: int, j: int) -> tuple[float, float]:
            return self.xa + i * hx, self.ya + j * hy

        K = [[0.0] * n_total for _ in range(n_total)]
        F = [0.0] * n_total

        for i in range(nx):
            for j in range(ny):
                x00, y00 = ncoords(i, j)
                x10, y10 = ncoords(i + 1, j)
                x01, y01 = ncoords(i, j + 1)
                x11, y11 = ncoords(i + 1, j + 1)

                self._add_element(
                    K,
                    F,
                    [(x00, y00), (x10, y10), (x01, y01)],
                    [nidx(i, j), nidx(i + 1, j), nidx(i, j + 1)],
                )
                self._add_element(
                    K,
                    F,
                    [(x10, y10), (x11, y11), (x01, y01)],
                    [nidx(i + 1, j), nidx(i + 1, j + 1), nidx(i, j + 1)],
                )

        boundary: dict[int, float] = {}
        for i in range(nx + 1):
            for j in range(ny + 1):
                if i == 0 or i == nx or j == 0 or j == ny:
                    k = nidx(i, j)
                    xk, yk = ncoords(i, j)
                    boundary[k] = self.g(xk, yk)

        for k, gk in boundary.items():
            for m in range(n_total):
                if m not in boundary:
                    F[m] -= K[m][k] * gk

        for k, gk in boundary.items():
            for m in range(n_total):
                K[k][m] = 0.0
                K[m][k] = 0.0
            K[k][k] = 1.0
            F[k] = gk

        A_mat = Matrix(*[Vector(*K[r]) for r in range(n_total)])
        b_vec = Vector(*F)
        sol = LinearSystem(A_mat, b_vec).solve(method=LinearSolverMethod.GAUSS_ELIMINATION)
        u = sol.components

        xa_ = self.xa
        ya_ = self.ya

        def eval_u(x: float, y: float) -> float:
            ix = max(0, min(nx - 1, int((x - xa_) / hx)))
            jy = max(0, min(ny - 1, int((y - ya_) / hy)))
            lx = max(0.0, min(1.0, (x - (xa_ + ix * hx)) / hx))
            ly = max(0.0, min(1.0, (y - (ya_ + jy * hy)) / hy))
            u00 = u[nidx(ix, jy)]
            u10 = u[nidx(ix + 1, jy)]
            u01 = u[nidx(ix, jy + 1)]
            u11 = u[nidx(ix + 1, jy + 1)]
            return (1 - lx) * (1 - ly) * u00 + lx * (1 - ly) * u10 + (1 - lx) * ly * u01 + lx * ly * u11

        return BivariateFunction(eval_u)

    def _add_element(self, K: list[list[float]], F: list[float], verts: list[tuple[float, float]], idxs: list[int]) -> None:
        x1, y1 = verts[0]
        x2, y2 = verts[1]
        x3, y3 = verts[2]

        area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
        if area == 0:
            return

        b = [y2 - y3, y3 - y1, y1 - y2]
        c = [x3 - x2, x1 - x3, x2 - x1]

        for i in range(3):
            for j in range(3):
                K[idxs[i]][idxs[j]] += (b[i] * b[j] + c[i] * c[j]) / (4 * area)
            F[idxs[i]] += (area / 3) * self.f(verts[i][0], verts[i][1])
