from typing import Literal
import math

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
                TOLERANCE=1e-10, N=100, 
                return_iterations=False, early_stop: int=None):

        if method == "bisection":
            assert a is not None, "a must be defined"
            assert b is not None, "b must be defined"
            assert a < b, "a must be less than b"
            assert self(a) * self(b) < 0, "f(a) and f(b) must have opposite signs"

            sol, n = self.bisection(a, b, TOLERANCE, N, early_stop)
            if return_iterations:
                return sol, n
            return sol
        
        if method == "newton":
            assert p0 is not None, "p0 must be defined"

            sol, n = self.newton(p0, TOLERANCE, N, early_stop)
            if return_iterations:
                return sol, n
            return sol
        
        if method == "secant":
            assert p0 is not None, "p0 must be defined"
            assert p1 is not None, "p1 must be defined"

            sol, n = self.secant(p0, p1, TOLERANCE, N, early_stop)
            if return_iterations:
                return sol, n
            return sol
        
        if method == "regula_falsi":
            assert p0 is not None, "p0 must be defined"
            assert p1 is not None, "p1 must be defined"
            assert self(p0) * self(p1) < 0, "f(p0) and f(p1) must have opposite signs"

            sol, n = self.regula_falsi(p0, p1, TOLERANCE, N, early_stop)
            if return_iterations:
                return sol, n
            return sol
        
        raise ValueError("Invalid method.")

    def bisection(self, a: float, b: float, TOLERANCE=1e-10, N=100, early_stop: int=None):
        for i in range(N):
            p = (a + b) / 2
            if self(p) == 0 or abs(a - b) < TOLERANCE or (early_stop is not None and i >= early_stop):
                return p, i + 1
            if self(a) * self(p) > 0:
                a = p
            else:
                b = p
        return None, N

    def newton(self, p0: float, TOLERANCE=1e-10, N=100, early_stop: int=None):
        deriv = self.differentiate()

        try:
            for i in range(N):
                p = p0 - self(p0) / deriv(p0)
                if abs(p - p0) < TOLERANCE or (early_stop is not None and i >= early_stop):
                    return p, i + 1
                p0 = p
            return None, N
        except ZeroDivisionError or OverflowError:
            return None, i
    
    def secant(self, p0: float, p1: float, TOLERANCE=1e-10, N=100, early_stop: int=None):
        for i in range(N):
            p = p1 - self(p1) * (p1 - p0) / (self(p1) - self(p0))
            if abs(p - p1) < TOLERANCE or (early_stop is not None and i >= early_stop):
                return p, i + 1
            p0 = p1
            p1 = p
        return None, N
    
    def regula_falsi(self, p0: float, p1: float, TOLERANCE=1e-10, N=100, early_stop: int=None):
        for i in range(N):
            p = p1 - self(p1) * (p1 - p0) / (self(p1) - self(p0))
            if abs(p - p1) < TOLERANCE or (early_stop is not None and i >= early_stop):
                return p, i + 1
            if self(p0) * self(p) > 0:
                p0 = p1
            p1 = p
        return None, N

    def fixed_point(self, p0: float, TOLERANCE=1e-10, N=100):
        assert p0 is not None, "p0 must be defined"

        try:
            for i in range(N):
                p = self(p0)
                if abs(p - p0) < TOLERANCE:
                    return p
                p0 = p
            return None
        except OverflowError:
            return None

class Polynomial(Function):

    def __init__(self, *coefficients):
        """
            coefficients are in the form a_0, a_1, ... a_n
        """
        self.function = lambda x: sum(a * x ** i for i, a in enumerate(coefficients))

class Exponent(Function):

    def __init__(self, f: Function, base: float=math.e):
        self.function = lambda x: base ** f(x)

class Sin(Function):

    def __init__(self, f: Function):
        self.function = lambda x: math.sin(f(x))

class Cos(Function):
    
    def __init__(self, f: Function):
        self.function = lambda x: math.cos(f(x))

class Tan(Function):
        
    def __init__(self, f: Function):
        self.function = lambda x: math.tan(f(x))
