from function import Function, Polynomial, Cos, Sin, Exponent, Log
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
        Approximate the integral of f(x) = x^3 5x^2 + 1 on the interval [1, 5]
        by using composite rectangular method and composite midpoint method. 
        """
        f = Polynomial(1, 0, 5, 1)
        
        # (a) with five subintervals
        res1 = f.integrate(1, 5, method='rectangular', n=5)
        res2 = f.integrate(1, 5, method='midpoint', n=5)
        
        # (b) with ten subintervals
        res3 = f.integrate(1, 5, method='rectangular', n=10)
        res4 = f.integrate(1, 5, method='midpoint', n=10)
        
        # (c) Compute the true error in both the cases
        f.integral(Polynomial(0, 1, 0, 5/3, 1/4))
        true = f.integrate(1, 5)
        
        return {
            'rectangular': {
                'n=5': res1,
                'n=10': res3,
                'true error with n=5': abs(true - res1),
                'true error with n=10': abs(true - res3)
            },
            'midpoint': {
                'n=5': res2,
                'n=10': res4,
                'true error with n=5': abs(true - res2),
                'true error with n=10': abs(true - res4)
            }
        }

    @staticmethod
    def problem2():
        """Redo Problem 1 by using composite trapezoidal method and composite Simpson method."""
        f = Polynomial(1, 0, 5, 1)
        
        # (a) with five subintervals
        res1 = f.integrate(1, 5, method='trapezoidal', n=5)
        res2 = f.integrate(1, 5, method='simpson', n=5)
        
        # (b) with ten subintervals
        res3 = f.integrate(1, 5, method='trapezoidal', n=10)
        res4 = f.integrate(1, 5, method='simpson', n=10)
        
        # (c) Compute the true error in both the cases
        f.integral(Polynomial(0, 1, 0, 5/3, 1/4))
        true = f.integrate(1, 5)
        
        return {
            'trapezoidal': {
                'n=5': res1,
                'n=10': res3,
                'true error with n=5': abs(true - res1),
                'true error with n=10': abs(true - res3)
            },
            'simpson': {
                'n=5': res2,
                'n=10': res4,
                'true error with n=5': abs(true - res2),
                'true error with n=10': abs(true - res4)
            }
        }

    @staticmethod
    def problem3():
        """
        Evaluate the following integral by using one point Gauss quadrature and
        compute the true error.
        int_0^{pi/2} x sin(x) dx
        """
        f = Polynomial(0, 1) * Sin(Polynomial(0, 1))
        f.integral(Sin(Polynomial(0, 1)) - Polynomial(0, 1) * Cos(Polynomial(0, 1)))
        
        res = f.integrate(0, math.pi/2, method='gauss', n=1)
        true = f.integrate(0, math.pi/2)
        
        return res, abs(true - res)
    
    @staticmethod
    def problem4():
        """Redo Problem 3 by using two point Gauss quadrature formula."""
        f = Polynomial(0, 1) * Sin(Polynomial(0, 1))
        f.integral(Sin(Polynomial(0, 1)) - Polynomial(0, 1) * Cos(Polynomial(0, 1)))
        
        res = f.integrate(0, math.pi/2, method='gauss', n=2)
        true = f.integrate(0, math.pi/2)
        
        return res, abs(true - res)        

    @staticmethod
    def problem5():
        """
        Use composite simpson's rule with n = 4 and m = 2 to approximate
        int_1.4^2.0 int_1.0^1.5 ln(x + 2y) dy dx
        """
        m, n = 2, 4
        a, b, c, d = 1.4, 2, 1, 1.5
        
        h = (d-c) / m
        f = lambda r: Log(Polynomial(2 * r, 1))
        
        g = (h / 6) * (
                f(d) + f(c) + 
                2 * sum(f(c + i * h) for i in range(1, m)) + 
                4 * sum(f(c + i * h + h/2) for i in range(m))
            )
        res = g.integrate(a, b, method='simpson', n=n)
        
        return res

    @staticmethod
    def problem6():
        pass
