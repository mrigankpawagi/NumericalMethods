from typing import Literal, Callable
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

    def differentiate(self, func=None, h=1e-5, method:Literal['forward', 'backward', 'central']='forward'):
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
            if method == 'forward':
                return self.differentiate_forward(h)
            if method == 'backward':
                return self.differentiate_forward(-h)
            if method == 'central':
                return self.differentiate_central(h)
            
            raise ValueError("Invalid method.")
        
    def differentiate_forward(self, h):
       return Function(lambda x: (self(x + h) - self(x)) / h)
   
    def differentiate_central(self, h):
        return Function(lambda x: (self(x + h) - self(x - h)) / (2 * h))
    
    def multi_differentiate(self, n: int, h=1e-5, method:Literal['forward', 'backward', 'central']='forward'):
        """
        Returns the nth derivative of the function.
        """
        if n == 0:
            return self
        return self.differentiate(h=h, method=method).multi_differentiate(n - 1, h, method)

    def integral(self, func=None, h=1e-5):
        """
        Sets or returns the integral of the function.
        If func is None, returns the integral.
        If func is a Function, sets the integral to func.
        If func is a lambda, sets the integral to a Function with the lambda.
        """
        if isinstance(func, Function):
            self.integral = func
        elif callable(func):
            self.integral = Function(func)
        else:
            raise NotImplementedError("Not implemented yet.")
        
    def integrate(self, a: float, b: float, method: Literal['rectangular', 'midpoint', 'trapezoidal', 'simpson', 'gauss']=None, n: int=None):
        """
        Definite integral of the function from a to b.
        """
        if method == 'rectangular':
            return self.integrate_rectangular(a, b, n)
        if method == 'midpoint':
            return self.integrate_midpoint(a, b, n)
        if method == 'trapezoidal':
            return self.integrate_trapezoidal(a, b, n)
        if method == 'simpson':
            return self.integrate_simpson(a, b, n)
        if method == 'gauss':
            return self.integrate_gauss(a, b, n)
        
        if hasattr(self, 'integral'):
            return self.integral(b) - self.integral(a)

        raise ValueError("Invalid method.")
    
    def integrate_rectangular(self, a: float, b: float, n: int = None):
        if not n: 
            return (b - a) * self(a)
        
        h = (b - a) / n
        return h * sum(self(a + i * h) for i in range(n))
    
    def integrate_midpoint(self, a: float, b: float, n: int = None):
        if not n: 
            return (b - a) * self((a + b) / 2)
        
        h = (b - a) / n
        return h * sum(self(a + i * h + h / 2) for i in range(n))
        
    def integrate_trapezoidal(self, a: float, b: float, n: int = None):
        if not n: 
            return (b - a) * (self(a) + self(b)) / 2
        
        h = (b - a) / n
        return h * (self(a) + 2 * sum(self(a + i * h) for i in range(1, n)) + self(b)) / 2
        
    def integrate_simpson(self, a: float, b: float, n: int = None):
        if not n: 
            return (b - a) * (self(a) + 4 * self((a + b) / 2) + self(b)) / 6
        
        h = (b - a) / n
        return h * (self(a) + 4 * sum(self(a + i * h + h / 2) for i in range(n)) + 2 * sum(self(a + i * h) for i in range(1, n)) + self(b)) / 6
    
    def integrate_gauss(self, a: float, b: float, n: int = None):
        t = Polynomial((a+b)/2, (b-a)/2)
        g = ((b-a)/2) * self(t) 
        if n == 1:
            return 2 * g(0)
        if n == 2:
            return g(-1/math.sqrt(3)) + g(1/math.sqrt(3))

        raise NotImplementedError("Not implemented except for n=1 and n=2.")
                        

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
        else:
            plt.show()


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
        
class MultiVariableFunction:

    def __init__(self, function):
        self.function = function
        
    def __call__(self, *args):
        unwrapped_args = []
        for arg in args:
            if isinstance(arg, Vector):
                unwrapped_args += arg.components
            else:
                unwrapped_args.append(arg)
        args = tuple(unwrapped_args)
        
        if all(arg is not None for arg in args):
            return self.function(*args)
        if all(arg is None for arg in args):
            raise ValueError("At least one argument must be defined.")

        def f(*z):
            arguments = []
            original_args = list(args)
            passed_args = list(z)
            for arg in original_args:
                if arg is None:
                    arguments.append(passed_args.pop(0))
                else:
                    arguments.append(arg)
            return self(*arguments)
        
        num_none = sum(arg is None for arg in args)
        if num_none == 1:
            return Function(f)
        return MultiVariableFunction(f)
        
class BivariateFunction(MultiVariableFunction):
    pass

class Vector:
    
    def __init__(self, *components):
        self.components = list(components)
    
    def __add__(self, other):
        return Vector(*[self.components[i] + other.components[i] for i in range(len(self.components))])
    
    def __sub__(self, other):
        return Vector(*[self.components[i] - other.components[i] for i in range(len(self.components))])
    
    def __rmul__(self, other):
        return Vector(*[other * self.components[i] for i in range(len(self.components))])
    
    def __call__(self, *args):
        return Vector(*[self.components[i](*args) for i in range(len(self.components))])
    
    def __getitem__(self, i):
        return self.components[i]
    
    def __setitem__(self, i, value):
        self.components[i] = value
    
    def __len__(self) -> int:
        return len(self.components)
    
    def __iter__(self):
        return iter(self.components)
    
    def __str__(self):
        return "<" + ", ".join(str(component) for component in self.components) + ">"
    
class Matrix:
        
    def __init__(self, *rows: list[Vector]):
        self.rows = list(rows)
    
    def __len__(self) -> int:
        return len(self.rows)
    
    def __getitem__(self, i):
        return self.rows[i]
    
    def __setitem__(self, i, value):
        self.rows[i] = value
    
    def __iter__(self):
        return iter(self.rows)
    
    def __add__(self, other):
        return Matrix(*[self.rows[i] + other.rows[i] for i in range(len(self))])
    
    def __sub__(self, other):
        return Matrix(*[self.rows[i] - other.rows[i] for i in range(len(self))])
    
    def __rmul__(self, other):
        return Matrix(*[other * self.rows[i] for i in range(len(self))])

class OrdinaryDifferentialEquation:

    def __init__(self):
        pass
    
class LinearODE(OrdinaryDifferentialEquation):

    def __init__(self):
        pass
    
class FirstOrderLinearODE(LinearODE):
    """
    y'(x) = f(x, y(x))
    These are initial value problems.
    """

    def __init__(self, f: BivariateFunction, a: float, b: float, y0: float):
        """
        f is a function of x and y(x)
        """
        self.f = f
        self.a = a
        self.b = b
        self.y0 = y0
        
    def solve(self, h: float = 0.1, method: Literal["euler", "runge-kutta", "taylor", "trapezoidal", "adam-bashforth", "adam-moulton", "predictor-corrector"]='euler', n: int = 1, step: int = 2, points: list[float]=[]):
        if method == "euler":
            return self.solve_taylor(h, 1)
        if method == "runge-kutta":
            return self.solve_runge_kutta(h, n)
        if method == "taylor":
            return self.solve_taylor(h, n)
        if method == "trapezoidal":
            return self.solve_trapezoidal(h)
        if method == "adam-bashforth":
            return self.solve_adam_bashforth(h, step, points)
        if method == "adam-moulton":
            return self.solve_adam_moulton(h, step, points)
        if method == "predictor-corrector":
            return self.solve_predictor_corrector(h)
        raise ValueError("Invalid method.")

    def solve_runge_kutta(self, h: float, n: int) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)
        if n == 1:
            return self.solve(h, method='euler')
        elif n == 2:
            for i in range(N):
                xi = self.a + i * h
                w.append(w[i] + (h/2) * self.f(xi, w[i]) + (h/2) * self.f(xi + h, w[i] + h * self.f(xi, w[i])))
        elif n == 3:
            for i in range(N):
                xi = self.a + i * h
                k1 = self.f(xi, w[i])
                k2 = self.f(xi + (h/3), w[i] + (h/3) * k1)
                k3 = self.f(xi + (2/3) * h, w[i] + (2/3) * h * k2)
                w.append(w[i] + (h/4) * (k1 + 3 * k3))
        elif n == 4:
            for i in range(N):
                xi = self.a + i * h
                k1 = h * self.f(xi, w[i])
                k2 = h * self.f(xi + h/2, w[i] + 0.5 * k1)
                k3 = h * self.f(xi + h/2, w[i] + 0.5 * k2)
                k4 = h * self.f(xi + h, w[i] + k3)
                w.append(w[i] + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4))
        else:
            raise NotImplementedError("Not implemented for n > 4 yet.")
        
        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except:
            return w # return list of values if interpolation fails (like when w is a list of Vectors)

    def solve_taylor(self, h: float, n: int) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)
        if n == 1:
            for i in range(N):
                xi = self.a + i * h
                w.append(w[i] + h * self.f(xi, w[i]))
        elif n == 2:
            for i in range(N):
                xi = self.a + i * h
                # g = f'
                g = self.f(None, w[i]).differentiate()(xi) + self.f(xi, w[i]) * self.f(xi, None).differentiate()(w[i])
                w.append(w[i] + h * self.f(xi, w[i]) + (h ** 2) * g / 2)
        else:
            raise NotImplementedError("Not implemented for n > 2 yet.")

        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except:
            return w # return list of values if interpolation fails (like when w is a list of Vectors)
    
    def solve_trapezoidal(self, h: float) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)
        for i in range(N):
            xi = self.a + i * h
            g = Function(lambda x: w[i] + (h/2) * (self.f(xi, w[i]) + self.f(xi + h, x)))
            w.append(g.fixed_point(w[i]))
            
        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except:
            return w # return list of values if interpolation fails (like when w is a list of Vectors)
    
    def solve_adam_bashforth(self, h: float, step: int, points: list[float]) -> Polynomial:
        w = [self.y0] + points
        N = int((self.b - self.a) / h)
        if step == 2:
            for i in range(1, N):
                xi = self.a + i * h
                w.append(w[i] + (h/2) * (3 * self.f(xi, w[i]) - self.f(xi - h, w[i-1])))
        elif step == 3:
            for i in range(2, N):
                xi = self.a + i * h
                w.append(w[i] + (h/12) * (23 * self.f(xi, w[i]) - 16 * self.f(xi - h, w[i-1]) + 5 * self.f(xi - 2 * h, w[i-2])))
        elif step == 4:
            for i in range(3, N):
                xi = self.a + i * h
                w.append(w[i] + (h/24) * (55 * self.f(xi, w[i]) - 59 * self.f(xi - h, w[i-1]) + 37 * self.f(xi - 2 * h, w[i-2]) - 9 * self.f(xi - 3 * h, w[i-3])))
        else:
            if step > 1:
                raise NotImplementedError("Not implemented for step > 4 yet.")
            else:
                raise ValueError("Step must be greater than 1.")
        
        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except:
            return w # return list of values if interpolation fails (like when w is a list of Vectors)
    
    def solve_adam_moulton(self, h: float, step: int, points: list[float]) -> Polynomial:
        w = [self.y0] + points
        N = int((self.b - self.a) / h)
        if step == 2:
            for i in range(1, N):
                xi = self.a + i * h
                g = Function(lambda x: w[i] + (h/12) * (5 * self.f(xi + h, x) + 8 * self.f(xi, w[i]) - self.f(xi - h, w[i-1])))
                w.append(g.fixed_point(w[i]))
        elif step == 3:
            for i in range(2, N):
                xi = self.a + i * h
                g = Function(lambda x: w[i] + (h/24) * (9 * self.f(xi + h, x) + 19 * self.f(xi, w[i]) - 5 * self.f(xi - h, w[i-1]) + self.f(xi - 2 * h, w[i-2])))
                w.append(g.fixed_point(w[i]))
        elif step == 4:
            for i in range(3, N):
                xi = self.a + i * h
                g = Function(lambda x: w[i] + (h/720) * (251 * self.f(xi + h, x) + 646 * self.f(xi, w[i]) - 264 * self.f(xi - h, w[i-1]) + 106 * self.f(xi - 2 * h, w[i-2]) - 19 * self.f(xi - 3 * h, w[i-3])))
                w.append(g.fixed_point(w[i]))
        else:
            if step > 1:
                raise NotImplementedError("Not implemented for step > 4 yet.")
            else:
                raise ValueError("Step must be greater than 1.")
        
        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except:
            return w # return list of values if interpolation fails (like when w is a list of Vectors)
    
    def solve_predictor_corrector(self, h: float) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)
        alphas = [self.a + h * i for i in range(1, 4)]

        # determine starting values with runge-kutta order-4
        g = self.f
        ivp = FirstOrderLinearODE(g, self.a, self.a + 3 * h, self.y0)
        rk_sol = ivp.solve_runge_kutta(h, 4)
        for i in range(3):
            w.append(rk_sol(alphas[i]))
            
        for i in range(3, N):
            xi = self.a + i * h
            
            # Predictor: Adams-Bashforth order-4
            prediction = w[i] + (h/24) * (55 * self.f(xi, w[i]) - 59 * self.f(xi - h, w[i-1]) + 37 * self.f(xi - 2 * h, w[i-2]) - 9 * self.f(xi - 3 * h, w[i-3]))
            
            # Corrector: Adams-Moulton order-3
            correction = w[i] + (h/24) * (9 * self.f(xi + h, prediction) + 19 * self.f(xi, w[i]) - 5 * self.f(xi - h, w[i-1]) + self.f(xi - 2 * h, w[i-2]))
            
            w.append(correction)
            
        try:
            return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
        except:
            return w # return list of values if interpolation fails (like when w is a list of Vectors)      

class SecondOrderLinearODE_BVP(LinearODE):
    """
    y''(x) = p(x)y'(x) + q(x)y(x) + r(x)
    These are boundary value problems.
    """

    def __init__(self, p: Function, q: Function, r: Function, a: float, b: float, y0: float, y1: float):
        self.p = p
        self.q = q
        self.r = r
        self.a = a
        self.b = b
        self.y0 = y0 # y(a)
        self.y1 = y1 # y(b)
        
    def solve(self, h: float = 0.1, method: Literal["shooting", "finite_difference"]="shooting"):
        if method == "shooting":
            return self.solve_shooting(h)
        if method == "finite_difference":
            return self.solve_finite_difference(h)
        raise ValueError("Invalid method.")
    
    def solve_shooting(self, h: float) -> Polynomial:
        IVP1 = SecondOrderODE_IVP(
            MultiVariableFunction(lambda t, u1, u2: self.p(t) * u2 + self.q(t) * u1 + self.r(t)),
            self.a, self.b, self.y0, 0
        )
        IVP2 = SecondOrderODE_IVP(
            MultiVariableFunction(lambda t, u1, u2: self.p(t) * u2 + self.q(t) * u1),
            self.a, self.b, 0, 1
        )
        
        sol1 = IVP1.solve(h)
        sol2 = IVP2.solve(h)
        
        c = (self.y1 - sol1(self.b)) / sol2(self.b)
        
        return sol1 + c * sol2
    
    def solve_finite_difference(self, h: float) -> Polynomial:
        N = int((self.b - self.a) / h) - 1
        A = Matrix(*[Vector(*[0 for _ in range(N)]) for _ in range(N)])
        b = Vector(*[0 for _ in range(N)])
        
        # First row
        A[0][0] = -(2 + (h ** 2) * self.q(self.a + h))
        A[0][1] = 1 - (h / 2) * self.p(self.a + h)
        b[0] = (h ** 2) * self.r(self.a + h) - (1 + (h / 2) * self.p(self.a + h)) * self.y0
        
        # Middle rows
        for i in range(1, N - 1):
            xi = self.a + (i+1) * h
            A[i][i-1] = 1 + (h / 2) * self.p(xi)
            A[i][i] = -(2 + (h ** 2) * self.q(xi))
            A[i][i+1] = 1 - (h / 2) * self.p(xi)
            b[i] = (h ** 2) * self.r(xi)
        
        # Last row
        A[N-1][N-2] = 1 + (h / 2) * self.p(self.b - h)
        A[N-1][N-1] = -(2 + (h ** 2) * self.q(self.b - h))
        b[N-1] = (h ** 2) * self.r(self.b - h) - (1 - (h / 2) * self.p(self.b - h)) * self.y1
        
        # Solve system of equations
        sol = LinearSystem(A, b).solve()
        
        w = [self.y0] + sol.components + [self.y1]
        
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 2)])
        

class SecondOrderODE_IVP(OrdinaryDifferentialEquation):
    """
    y''(x) = f(x, y(x), y'(x))
    These are initial value problems.
    """
    
    def __init__(self, f: MultiVariableFunction, a: float, b: float, y0: float, y1: float):
        """
        f is a function of x, y(x), and y'(x)
        """
        self.f = f
        self.a = a
        self.b = b
        self.y0 = y0 # y(a)
        self.y1 = y1 # y'(a)
        
    def solve(self, h: float = 0.1, method: Literal["euler", "runge-kutta", "taylor", "trapezoidal", "adam-bashforth", "adam-moulton", "predictor-corrector"]='euler', n: int = 1, step: int = 2, points: list[float]=[]):
        U0 = Vector(self.y0, self.y1)
        F =  Vector(
            MultiVariableFunction(lambda t, u1, u2: u2),
            MultiVariableFunction(lambda t, u1, u2: self.f(t, u1, u2))
        )
        IVP = FirstOrderLinearODE(F, self.a, self.b, U0)
        sol = IVP.solve(h, method, n, step, points)
        
        w = [x[0] for x in sol]
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(len(w))])

class SecondOrderODE_BVP(OrdinaryDifferentialEquation):
    """
    y''(x) = f(x, y(x), y'(x))
    These are boundary value problems.
    """
    
    def __init__(self, f: MultiVariableFunction, a: float, b: float, y0: float, y1: float):
        """
        f is a function of x, y(x), and y'(x)
        """
        self.f = f
        self.a = a
        self.b = b
        self.y0 = y0 # y(a)
        self.y1 = y1 # y(b)
        
    def solve(self, h: float = 0.1, method: Literal["shooting_newton"]="shooting_newton", M: int = 100, TOL: float = 1e-5, initial_approximation=None):
        if method == "shooting_newton":
            return self.solve_shooting_newton(h, M, TOL, initial_approximation)
        raise ValueError("Invalid method.")
    
    def solve_shooting_newton(self, h: float, M, TOL, initial_approximation) -> Polynomial:
        t = 1 if initial_approximation is None else initial_approximation # initial guess for y'(a)
        i = 0
        
        while i < M:
            IVP1 = SecondOrderODE_IVP(
                MultiVariableFunction(lambda t, u1, u2: self.f(t, u1, u2)),
                self.a, self.b, self.y0, t
            )
            y = IVP1.solve(h)
            
            p = Function(lambda x: self.f(x, None, y.differentiate()(x)).differentiate()(y(x)))
            q = Function(lambda x: self.f(x, y(x), None).differentiate()(y.differentiate()(x)))
            r = Function(lambda x: 0)
            IVP2 = SecondOrderLinearODE_BVP(p, q, r, self.a, self.b, 0, 1)
            z = IVP2.solve(h)
            
            t0 = t - (y(self.b) - self.y1) / z(self.b)
            if abs(t0 - t) < TOL:
                return y
            
            t = t0            
            i += 1
        
        return None

class LinearSystem:
    """
    A system of linear equations.
    """
    
    def __init__(self, A: Matrix, b: Vector):
        self.A = A
        self.b = b
        self.N = len(A)
        
        assert len(b) == self.N, "A and b must have the same number of rows."
        
    def solve(self, method: Literal["gauss_elimination", "gauss_jacobi", "gauss_seidel"]='gauss_elimination', TOL: float = 1e-5, initial_approximation: Vector = None, MAX_ITERATIONS: int = 100):
        if method == "gauss_elimination":
            return self.solve_gauss_elimination()
        if method == "gauss_jacobi":
            return self.solve_gauss_jacobi(TOL, initial_approximation, MAX_ITERATIONS)
        if method == "gauss_seidel":
            return self.solve_gauss_seidel(TOL, initial_approximation, MAX_ITERATIONS)
        raise ValueError("Invalid method.")
    
    def solve_gauss_elimination(self) -> Vector:
        for i in range(self.N - 1):
            # find pivot row
            p = None
            for j in range(i, self.N):
                if abs(self.A[j][i]) != 0:
                    p = j
                    break
            if p is None:
                raise ValueError("No unique solution exists.")
            
            if p != i:
                # swap rows
                self.A[i], self.A[p] = self.A[p], self.A[i]
                self.b[i], self.b[p] = self.b[p], self.b[i]
            
            for j in range(i + 1, self.N):
                m = self.A[j][i] / self.A[i][i]
                self.A[j] = self.A[j] - m * self.A[i]
                self.b[j] = self.b[j] - m * self.b[i]
            
        if abs(self.A[self.N - 1][self.N - 1]) == 0:
            raise ValueError("No unique solution exists.")
        
        x = [0] * self.N
        x[self.N - 1] = self.b[self.N - 1] / self.A[self.N - 1][self.N - 1]
        for i in range(self.N - 2, -1, -1):
            x[i] = (self.b[i] - sum(self.A[i][j] * x[j] for j in range(i + 1, self.N))) / self.A[i][i]
        
        return Vector(*x)
    
    def solve_gauss_jacobi(self, TOL: float, initial_approximation: Vector, MAX_ITERATIONS) -> Vector:
        assert initial_approximation is not None, "Initial approximation must be defined."
        
        x0 = initial_approximation
        k = 0
        
        while k < MAX_ITERATIONS:
            x1 = [0] * self.N
            for i in range(self.N):
                x1[i] = (self.b[i] - sum(self.A[i][j] * x0[j] for j in range(self.N) if j != i)) / self.A[i][i]
            
            if max(abs(x1[i] - x0[i]) for i in range(self.N)) / max(abs(x1[i]) for i in range(self.N)) < TOL:
                return Vector(*x1)
            
            x0 = x1
            k += 1
            
        return None
    
    def solve_gauss_seidel(self, TOL: float, initial_approximation: Vector, MAX_ITERATIONS) -> Vector:
        assert initial_approximation is not None, "Initial approximation must be defined."
        
        x0 = initial_approximation
        k = 0
        
        while k < MAX_ITERATIONS:
            x1 = [0] * self.N
            for i in range(self.N):
                x1[i] = (self.b[i] - sum(self.A[i][j] * x1[j] for j in range(i)) - sum(self.A[i][j] * x0[j] for j in range(i+1, self.N))) / self.A[i][i]
                
            if max(abs(x1[i] - x0[i]) for i in range(self.N)) / max(abs(x1[i]) for i in range(self.N)) < TOL:
                return Vector(*x1)
            
            x0 = x1
            k += 1
            
        return None
