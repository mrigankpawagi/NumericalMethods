from __future__ import annotations

from .vector import Vector


class Matrix:
    """A rectangular array of scalars stored as a sequence of row :class:`Vector` objects."""

    def __init__(self, *rows: Vector) -> None:
        self.rows: list[Vector] = list(rows)

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other: Matrix) -> Matrix:
        return Matrix(*[self.rows[i] + other.rows[i] for i in range(len(self))])

    def __sub__(self, other: Matrix) -> Matrix:
        return Matrix(*[self.rows[i] - other.rows[i] for i in range(len(self))])

    def __rmul__(self, scalar: float) -> Matrix:
        return Matrix(*[scalar * self.rows[i] for i in range(len(self))])

    # ------------------------------------------------------------------
    # Container interface
    # ------------------------------------------------------------------

    def __getitem__(self, i: int) -> Vector:
        return self.rows[i]

    def __setitem__(self, i: int, value: Vector) -> None:
        self.rows[i] = value

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __repr__(self) -> str:
        return f"Matrix({', '.join(repr(r) for r in self.rows)})"
