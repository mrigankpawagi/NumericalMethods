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
        
class BivariateFunction:

    def __init__(self, function):
        self.function = function
        
    def __call__(self, x=None, y=None):
        if x is not None and y is not None:
            return self.function(x, y)
        if x is not None:
            return Function(lambda z: self(x, z))
        if y is not None:
            return Function(lambda z: self(z, y))
        raise ValueError("x or y must be defined.")

class OrdinaryDifferentialEquation:

    def __init__(self):
        pass
    
class LinearODE(OrdinaryDifferentialEquation):

    def __init__(self):
        pass
    
class FirstOrderLinearODE(LinearODE):
    """
    y'(x) = f(x, y(x))
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
        
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
    
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

        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
    
    def solve_trapezoidal(self, h: float) -> Polynomial:
        w = [self.y0]
        N = int((self.b - self.a) / h)
        for i in range(N):
            xi = self.a + i * h
            g = Function(lambda x: w[i] + (h/2) * (self.f(xi, w[i]) + self.f(xi + h, x)))
            w.append(g.fixed_point(w[i]))
            
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
    
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
        
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
    
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
        
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])
    
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
            
        return Polynomial.interpolate([(self.a + i * h, w[i]) for i in range(N + 1)])           
