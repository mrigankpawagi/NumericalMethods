from function import Polynomial, Log, BivariateFunction, FirstOrderLinearODE, Exponent
from util import Util
import math

class Solution:

    @staticmethod
    def solve(n):
        return getattr(Solution, f'problem{n}')()

    @staticmethod
    def solve_all(n):
        ans = [Solution.solve(i) for i in range(1, n + 1)]
        for i in range(len(ans)):
            print(f"Problem {i+1}: {ans[i]}")

    @staticmethod
    def problem1():
        """
        Use the forward difference formula to approximate the derivative of f(x) = ln x
        at x_0 = 1.8 using h = 0.1, 0.05, 0.01. Determine the bounds for the approximation errors.
        """
        f = Log(Polynomial(0, 1))
        a = 1.8
        max_err_func = lambda h: h / (2 * (a**2))
        
        # h = 0.1
        res1 = f.differentiate(h=0.1, method='forward')(a)
        b1 = max_err_func(0.1) # error bound
        
        # h = 0.05
        res2 = f.differentiate(h=0.05, method='forward')(a)
        b2 = max_err_func(0.05) # error bound
        
        # h = 0.01
        res3 = f.differentiate(h=0.01, method='forward')(a)
        b3 = max_err_func(0.01) # error bound
        
        return {
            'h=0.1': {
                'result': res1,
                'error_bound': b1
            },
            'h=0.05': {
                'result': res2,
                'error_bound': b2
            },
            'h=0.001': {
                'result': res3,
                'error_bound': b3
            }
        }

    @staticmethod
    def problem2():
        """
        Redo problem 1 with the backward difference formula and the central difference formula.
        """
        f = Log(Polynomial(0, 1))
        a = 1.8

        # backward difference
        max_err_func = lambda h: h / (2 * ((a - h) ** 2))
        
        # h = 0.1
        res1 = f.differentiate(h=0.1, method='backward')(a)
        
        # h = 0.05
        res2 = f.differentiate(h=0.05, method='backward')(a)
        
        # h = 0.01
        res3 = f.differentiate(h=0.01, method='backward')(a)
        
        backward = {
            'h=0.1': {
                'result': res1,
                'error_bound': max_err_func(0.1)
            },
            'h=0.05': {
                'result': res2,
                'error_bound': max_err_func(0.05)
            },
            'h=0.001': {
                'result': res3,
                'error_bound': max_err_func(0.01)
            }
        }
        
        # central difference
        max_err_func = lambda h: (h**2) / (3 * ((a-h)**3))
        
        # h = 0.1
        res1 = f.differentiate(h=0.1, method='central')(a)
        
        # h = 0.05
        res2 = f.differentiate(h=0.05, method='central')(a)
        
        # h = 0.01
        res3 = f.differentiate(h=0.01, method='central')(a)
        
        central = {
            'h=0.1': {
                'result': res1,
                'error_bound': max_err_func(0.1)
            },
            'h=0.05': {
                'result': res2,
                'error_bound': max_err_func(0.05)
            },
            'h=0.001': {
                'result': res3,
                'error_bound': max_err_func(0.01)
            }
        }
        
        return {
            'backward': backward,
            'central': central
        }

    @staticmethod
    def problem3():
        """
        Consider the IVP y' = y ln y / x, y(2) = e. 
        Use Euler's method with h = 0.1 to obtain the approximation to y(3).
        """
        f = BivariateFunction(lambda x, y: Polynomial(0, 1)(y) * Log(Polynomial(0, 1))(y) / Polynomial(0, 1)(x))
        a = 2
        b = 3
        IVP = FirstOrderLinearODE(f, a, b, math.e)
        sol = IVP.solve(h=0.1, method='euler')
        
        return sol(b)
    
    @staticmethod
    def problem4():
        """
        Consider the IVP y' = y - x, y(0) = 1/2.
        Use Euler's method with h = 0.1 and h = 0.05 to obtain the approximation to y(1).
        Given that the exact solution to the IVP is y(x) = x + 1 - 1/2 e^x, compare the errors in the two approximations.
        """
        a = 0
        b = 1
        f = BivariateFunction(lambda x, y: y - x)
        gt = (Polynomial(1, 1) - (0.5 * Exponent(Polynomial(0, 1))))(b)
        IVP = FirstOrderLinearODE(f, a, b, 0.5)
        
        # h = 0.1
        sol1 = IVP.solve(h=0.1, method='euler')(b)
        err1 = abs(gt - sol1)
        
        h = 0.05
        sol2 = IVP.solve(h=0.05, method='euler')(b)
        err2 = abs(gt - sol2)
        
        return {
            'h=0.1': {
                'result': sol1,
                'error': err1
            },
            'h=0.05': {
                'result': sol2,
                'error': err2
            }
        }        

    @staticmethod
    def problem5():
        """
        Consider the IVP y' = 2xy^2, y(0) = 0.5.
        Use Euler's method with h = 0.1 to obtain the approximation to y(1).
        """
        a = 0
        b = 1
        f = BivariateFunction(lambda x, y: 2 * x * (y**2))
        IVP = FirstOrderLinearODE(f, a, b, 0.5)
        sol = IVP.solve(h=0.03, method='euler')
        
        return sol(b)
