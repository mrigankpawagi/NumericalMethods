from typing import Literal

class Function:

    def __init__(self, function):
        self.function = function

    def __call__(self, x):
        return self.function(x)

    def __truediv__(f, g):
        return Function(lambda x: f(x) / g(x))
    
    def __mul__(f, g):
        return Function(lambda x: f(x) * g(x))
    
    def __add__(f, g):
        return Function(lambda x: f(x) + g(x))
    
    def __sub__(f, g):
        return Function(lambda x: f(x) - g(x))

    def differentiate(self, func=None, h=10e-5):
        """
        Sets or returns the derivative of the function.
        If func is None, returns the derivative.
        If func is a Function, sets the derivative to func.
        If func is a lambda, sets the derivative to a Function with the lambda.
        """
        if isinstance(func, Function):
            self.derivative = func
        elif callable(func):
            self.derivative = Function(func)
        else:
            return Function(lambda x: (self(x + h) - self(x)) / h)

    def root(self, method: Literal["bisection", "newton", "secant", "regula_falsi"],
                a: float = None, b: float = None, 
                p0: float = None, p1: float = None,
                TOLERANCE=1e-10, N=100):

        if method == "bisection":
            assert a is not None, "a must be defined"
            assert b is not None, "b must be defined"
            assert a < b, "a must be less than b"
            assert self(a) * self(b) < 0, "f(a) and f(b) must have opposite signs"

            return self.bisection(a, b, TOLERANCE, N)
        
        if method == "newton":
            assert p0 is not None, "p0 must be defined"

            return self.newton(p0, TOLERANCE, N)
        
        if method == "secant":
            assert p0 is not None, "p0 must be defined"
            assert p1 is not None, "p1 must be defined"

            return self.secant(p0, p1, TOLERANCE, N)
        
        if method == "regula_falsi":
            assert p0 is not None, "p0 must be defined"
            assert p1 is not None, "p1 must be defined"
            assert self(p0) * self(p1) < 0, "f(p0) and f(p1) must have opposite signs"

            return self.regula_falsi(p0, p1, TOLERANCE, N)

    def bisection(self, a: float, b: float, TOLERANCE=1e-10, N=100):
        for _ in range(N):
            p = (a + b) / 2
            if self(p) == 0 or abs(a - b) < TOLERANCE:
                return p
            if self(a) * self(p) > 0:
                a = p
            else:
                b = p
        return None

    def newton(self, p0: float, TOLERANCE=1e-10, N=100):
        deriv = self.differentiate()

        for _ in range(N):
            p = p0 - self(p0) / deriv(p0)
            if abs(p - p0) < TOLERANCE:
                return p
            p0 = p
        return None
    
    def secant(self, p0: float, p1: float, TOLERANCE=1e-10, N=100):
        for _ in range(N):
            p = p1 - self(p1) * (p1 - p0) / (self(p1) - self(p0))
            if abs(p - p1) < TOLERANCE:
                return p
            p0 = p1
            p1 = p
        return None
    
    def regula_falsi(self, p0: float, p1: float, TOLERANCE=1e-10, N=100):
        for _ in range(N):
            p = p1 - self(p1) * (p1 - p0) / (self(p1) - self(p0))
            if abs(p - p1) < TOLERANCE:
                return p
            if self(p0) * self(p) > 0:
                p0 = p1
            p1 = p
        return None

    def fixed_point(self, p0: float, TOLERANCE=1e-10, N=100):
        assert p0 is not None, "p0 must be defined"

        for _ in range(N):
            p = self(p0)
            if abs(p - p0) < TOLERANCE:
                return p
            p0 = p
        return None

class Polynomial(Function):

    def __init__(self, *coefficients):
        """
            coefficients are in the form a_0, a_1, ... a_n
        """
        self.function = lambda x: sum(a * x ** i for i, a in enumerate(coefficients))
