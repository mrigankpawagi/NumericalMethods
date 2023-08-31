from typing import Literal
import math

class Function:

    def __init__(self, function):
        self.function = function

    def __call__(self, x):
        if callable(x):
            return Function(lambda p: self(x(p)))
        return self.function(x)

    def __truediv__(f, g):
        return Function(lambda x: f(x) / g(x))
    
    def __mul__(f, g):
        return Function(lambda x: f(x) * g(x))
    
    def __add__(f, g):
        return Function(lambda x: f(x) + g(x))
    
    def __sub__(f, g):
        return Function(lambda x: f(x) - g(x))
    
    def __pow__(f, g):
        return Function(lambda x: f(x) ** g(x))
    
    def __rtruediv__(f, g):
        return Function(lambda x: g / f(x))
    
    def __rmul__(f, g):
        return Function(lambda x: g * f(x))
    
    def __radd__(f, g):
        return Function(lambda x: g + f(x))
    
    def __rsub__(f, g):
        return Function(lambda x: g - f(x))
    
    def __neg__(f):
        return Function(lambda x: -f(x))

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

    def root(self, method: Literal["bisection", "newton", "secant", "regula_falsi", "modified_newton"],
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
        
        if method == "modified_newton":
            assert p0 is not None, "p0 must be defined"

            sol, n = self.modified_newton(p0, TOLERANCE, N, early_stop)
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
        
    def modified_newton(self, p0: float, TOLERANCE=1e-10, N=100, early_stop: int=None):
        deriv = self.differentiate()
        double_deriv = deriv.differentiate()

        try:
            for i in range(N):
                p = p0 - self(p0) * deriv(p0) / (deriv(p0) ** 2 - self(p0) * double_deriv(p0))
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
        
    def plot(self, min: float, max: float, N=1000, file: str="", clear: bool=False):
        import numpy as np
        import matplotlib.pyplot as plt

        # get N equally spaced points in [min, max]
        x = [min + (i/N) * + (max - min) for i in range(N)]
        y = [self(t) for t in x]

        if clear:
            plt.clf()
        plt.plot(x, y)
        plt.xlabel('x')
        plt.ylabel('y')
        if file:
            plt.savefig(file)


class Polynomial(Function):

    def __init__(self, *coefficients):
        """
            coefficients are in the form a_0, a_1, ... a_n
        """
        self.function = lambda x: sum(a * x ** i for i, a in enumerate(coefficients))
        self.coefficients = coefficients

    @staticmethod
    def interpolate(data: tuple, method: Literal["lagrange", "newton"]='newton', f: Function=None,
                    form: Literal["standard", "backward_diff", "forward_diff"]='standard'):
        """
            data is a list of (x, y) tuples.
            alternative: f is a Function that returns the y values and data is a list of x values.
        """
        if f is not None:
            data = [(x, f(x)) for x in data]

        if method == "lagrange":
            return Polynomial.interpolate_lagrange(data)
        if method == "newton":
            if form == "standard":
                return Polynomial.interpolate_newton(data)
            if form == "backward_diff":
                return Polynomial.interpolate_newton_backward_diff(data)
            if form == "forward_diff":
                return Polynomial.interpolate_newton_forward_diff(data)
        raise ValueError("Invalid method.")
    
    @staticmethod
    def interpolate_lagrange(data: tuple):
        """
            data is a tuple of (x, y) tuples
        """
        n = len(data)
        x = [data[i][0] for i in range(n)]
        y = [data[i][1] for i in range(n)]
        p = Polynomial(0)

        for i in range(n):
            L = Polynomial(1)
            for j in range(n):
                if i != j:
                    L *= Polynomial(-x[j], 1) / Polynomial(x[i] - x[j])
            p += y[i] * L

        return p
    
    @staticmethod
    def interpolate_newton(data: tuple):
        """
            data is a tuple of (x, y) tuples
        """
        n = len(data)
        x = [data[i][0] for i in range(n)]
        y = [data[i][1] for i in range(n)]

        def divided_difference(i, j):
            if i == j:
                return y[i]
            return (divided_difference(i + 1, j) - divided_difference(i, j - 1)) / (x[j] - x[i])
        
        def factor_product(roots: list):
            if not roots:
                return Polynomial(1)
            return Polynomial(-roots[0], 1) * factor_product(roots[1:])
        
        coefficients = [divided_difference(0, i) for i in range(n)]
        
        p = Polynomial(coefficients[0])
        for i in range(1, n):
            p = p + coefficients[i] * factor_product(x[:i])

        return p
    
    @staticmethod
    def interpolate_newton_backward_diff(data: tuple):
        """
            data is a tuple of (x, y) tuples
        """
        from util import Util
        data = sorted(data, key=lambda x: x[0])
        diffs = sorted([data[i+1][0] - data[i][0] for i in range(len(data) - 1)])
        for j in range(len(diffs) - 1):
            assert abs(diffs[j+1] - diffs[j]) < 1e-6, "x values must be equally spaced"

        h = abs(diffs[0])
        n = len(data) - 1
        
        p = Polynomial(data[n][1])
        for k in range(1, n+1):
            p += (((-1) ** k) * Util.downdelta(k, n, data)) * Util.choose(Polynomial(data[n][0] / h, - 1 / h), k)

        return p       


    @staticmethod
    def interpolate_newton_forward_diff(data: tuple):
        """
            data is a tuple of (x, y) tuples
        """
        from util import Util
        data = sorted(data, key=lambda x: x[0])
        diffs = sorted([data[i+1][0] - data[i][0] for i in range(len(data) - 1)])
        for j in range(len(diffs) - 1):
            assert abs(diffs[j+1] - diffs[j]) < 1e-6, "x values must be equally spaced"

        h = abs(diffs[0])
        n = len(data) - 1
        
        p = Polynomial(data[0][1])
        for k in range(1, n+1):
            p += Util.delta(k, 0, data) * Util.choose(Polynomial(-data[0][0] / h, 1 / h), k)

        return p       


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

class Log(Function):

    def __init__(self, f: Function, base: float=math.e):
        self.function = lambda x: math.log(f(x), base)
