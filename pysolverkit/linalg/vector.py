from __future__ import annotations


class Vector:
    """An ordered sequence of scalar components with element-wise arithmetic."""

    def __init__(self, *components: float) -> None:
        self.components: list[float] = list(components)

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------

    def __add__(self, other: Vector) -> Vector:
        return Vector(*[self.components[i] + other.components[i] for i in range(len(self))])

    def __sub__(self, other: Vector) -> Vector:
        return Vector(*[self.components[i] - other.components[i] for i in range(len(self))])

    def __rmul__(self, scalar: float) -> Vector:
        return Vector(*[scalar * self.components[i] for i in range(len(self))])

    def __call__(self, *args) -> Vector:
        return Vector(*[self.components[i](*args) for i in range(len(self))])

    # ------------------------------------------------------------------
    # Container interface
    # ------------------------------------------------------------------

    def __getitem__(self, i: int) -> float:
        return self.components[i]

    def __setitem__(self, i: int, value: float) -> None:
        self.components[i] = value

    def __len__(self) -> int:
        return len(self.components)

    def __iter__(self):
        return iter(self.components)

    def __str__(self) -> str:
        return "<" + ", ".join(str(c) for c in self.components) + ">"

    def __repr__(self) -> str:
        return f"Vector({', '.join(repr(c) for c in self.components)})"
