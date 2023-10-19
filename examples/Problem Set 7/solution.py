from function import Function, Polynomial, Cos, Sin, Exponent, Log, FirstOrderLinearODE, BivariateFunction
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
        """Use Taylor's series method of order 2 to approximate the solution for each
        of the following initial value problems.
        """
        def answer(f, a, b, y0, h):
            return FirstOrderLinearODE(f, a, b, y0).solve(h, method='taylor', n=2)(b)

        # (a) y' = y/x - (y/x)^2, y(1) = 1, 1 <= x <= 2, h = 0.1
        ans1 = answer(BivariateFunction(lambda x, y: y / x - (y / x)**2), 1, 2, 1, 0.1)
        
        # (b) y' = sinx + e^-x, y(0) = 0, 0 <= x <= 1, h = 0.5
        ans2 = answer(BivariateFunction(lambda x, y: Sin(Polynomial(0,1))(x) + Exponent(Polynomial(0, -1))(x)), 0, 1, 0, 0.5)
        
        # (c) y' = (y^2 + y)/x, y(1) = -2, 1 <= x <= 3, h = 0.5
        ans3 = answer(BivariateFunction(lambda x, y: (y**2 + y) / x), 1, 3, -2, 0.5)
        
        # (d) y' = -xy + 4x/y, y(0) = 1, 0 <= x <= 1, h = 0.25
        ans4 = answer(BivariateFunction(lambda x, y: -x*y + 4*x/y), 0, 1, 1, 0.25)
        
        return ans1, ans2, ans3, ans4

    @staticmethod
    def problem2():
        """
        Redo Problem 1 using the Runga Kutta method of order 2.
        """
        def answer(f, a, b, y0, h):
            return FirstOrderLinearODE(f, a, b, y0).solve(h, method='rk', n=2)(b)

        # (a) y' = y/x - (y/x)^2, y(1) = 1, 1 <= x <= 2, h = 0.1
        ans1 = answer(BivariateFunction(lambda x, y: y / x - (y / x)**2), 1, 2, 1, 0.1)
        
        # (b) y' = sinx + e^-x, y(0) = 0, 0 <= x <= 1, h = 0.5
        ans2 = answer(BivariateFunction(lambda x, y: Sin(Polynomial(0,1))(x) + Exponent(Polynomial(0, -1))(x)), 0, 1, 0, 0.5)
        
        # (c) y' = (y^2 + y)/x, y(1) = -2, 1 <= x <= 3, h = 0.5
        ans3 = answer(BivariateFunction(lambda x, y: (y**2 + y) / x), 1, 3, -2, 0.5)
        
        # (d) y' = -xy + 4x/y, y(0) = 1, 0 <= x <= 1, h = 0.25
        ans4 = answer(BivariateFunction(lambda x, y: -x*y + 4*x/y), 0, 1, 1, 0.25)
        
        return ans1, ans2, ans3, ans4

    @staticmethod
    def problem3():
        """
        Redo Problem 1 using the Trapezoidal method.
        """
        def answer(f, a, b, y0, h):
            return FirstOrderLinearODE(f, a, b, y0).solve(h, method='trapezoidal')(b)

        # (a) y' = y/x - (y/x)^2, y(1) = 1, 1 <= x <= 2, h = 0.1
        ans1 = answer(BivariateFunction(lambda x, y: y / x - (y / x)**2), 1, 2, 1, 0.1)
        
        # (b) y' = sinx + e^-x, y(0) = 0, 0 <= x <= 1, h = 0.5
        ans2 = answer(BivariateFunction(lambda x, y: Sin(Polynomial(0,1))(x) + Exponent(Polynomial(0, -1))(x)), 0, 1, 0, 0.5)
        
        # (c) y' = (y^2 + y)/x, y(1) = -2, 1 <= x <= 3, h = 0.5
        ans3 = answer(BivariateFunction(lambda x, y: (y**2 + y) / x), 1, 3, -2, 0.5)
        
        # (d) y' = -xy + 4x/y, y(0) = 1, 0 <= x <= 1, h = 0.25
        ans4 = answer(BivariateFunction(lambda x, y: -x*y + 4*x/y), 0, 1, 1, 0.25)
        
        return ans1, ans2, ans3, ans4
        
    
    @staticmethod
    def problem4():
        """
        Using the Taylor's series method of order 2, Runga Kutta method of order 2, and Trapezoidal method to approximate the solution of the following initial value problems and compare the
        results.
        """
        def answer(f, a, b, y0, h, gt):
            sol1 = FirstOrderLinearODE(f, a, b, y0).solve(h, method='taylor', n=2)(b)
            err1 = abs(gt(b) - sol1)
        
            sol2 = FirstOrderLinearODE(f, a, b, y0).solve(h, method='rk', n=2)(b)
            err2 = abs(gt(b) - sol2)
            
            sol3 = FirstOrderLinearODE(f, a, b, y0).solve(h, method='trapezoidal')(b)
            err3 = abs(gt(b) - sol3)
            
            return {
                'taylor': {
                    'result': sol1,
                    'error': err1
                },
                'rk': {
                    'result': sol2,
                    'error': err2
                },
                'trapezoidal': {
                    'result': sol3,
                    'error': err3
                }
            }
        
        # (a) y' = xe^(3x) - 2y, y(0) = 0, 0 <= x <= 1, h = 0.1
        # actual solution: y = (1/5)xe^(3x) - (1/25)e^(3x) + (1/25)e^(-2x)
        ans1 = answer(BivariateFunction(lambda x, y: x * Exponent(Polynomial(0, 3))(x) - 2 * y),
                      0, 1, 0, 0.1,
                       Polynomial(0, 1/5) * Exponent(Polynomial(0, 3)) - (1/25) * Exponent(Polynomial(0, 3)) + (1/25) * Exponent(Polynomial(0, -2)))
        
        # (b) y' = 1 + (x-y)^2, y(2) = 1, 2 <= x <= 3, h = 0.5
        # actual solution: y = x + 1/(1-x)
        ans2 = answer(BivariateFunction(lambda x, y: 1 + (x - y)**2),
                      2, 3, 1, 0.5,
                      Polynomial(0, 1) + (1 / (Polynomial(1, -1))))
        
        # (c) y' = 1 + y/x, y(1) = 1, 1 <= x <= 2, h = 0.25
        # actual solution: y = xlnx + 2x
        ans3 = answer(BivariateFunction(lambda x, y: 1 + y / x),
                      1, 2, 1, 0.25,
                      Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2))
        
        return {
            'a': ans1,
            'b': ans2,
            'c': ans3
        }
