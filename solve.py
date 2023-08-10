from typing import Callable, Literal
from function import Function

class Solve:

    @staticmethod
    def root(f: Callable, 
                method: Literal["bisection", "fixed_point"],
                a: float = None, b: float = None, p0: float = None, 
                TOLERANCE=1e-10, N=100):

        if method == "bisection":
            assert a is not None, "a must be defined"
            assert b is not None, "b must be defined"
            assert a < b, "a must be less than b"
            assert f(a) * f(b) < 0, "f(a) and f(b) must have opposite signs"

            return Solve.bisection(f, a, b, TOLERANCE, N)

        if method == "fixed_point":
            assert p0 is not None, "p0 must be defined"

            return Solve.fixed_point(f, p0, TOLERANCE, N)

    @staticmethod
    def bisection(f, a, b, TOLERANCE, N):
        for _ in range(N):
            p = (a + b) / 2
            if f(p) == 0 or abs(a - b) < TOLERANCE:
                return p
            if f(a) * f(p) > 0:
                a = p
            else:
                b = p
        return None
    
    def fixed_point(f, p0, TOLERANCE, N):
        g = Function.FixedPoint(f)

        for _ in range(N):
            p = g(p0)
            if abs(p - p0) < TOLERANCE:
                return p
            p0 = p
        return None
