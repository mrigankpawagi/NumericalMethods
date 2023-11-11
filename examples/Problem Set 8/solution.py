from function import Function, Polynomial, Cos, Sin, Exponent, Log, FirstOrderLinearODE, BivariateFunction, Vector, MultiVariableFunction
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
        Use two step Adam-Bashforth explicit method to approximate
        the solutions of the following initial value problems. Compute
        the value of the solution at the end point of the interval and
        find the error.
        """
        def answer(f, a, b, y0, h, gt):
            sol = FirstOrderLinearODE(f, a, b, y0).solve(h, method='adam-bashforth', step=2, points=[gt(a+h)])(b)
            return {
                'solution': sol,
                'error': abs(gt(b) - sol)
            }
        
        # (a) y' = t e^(3t) - 2y, y(0) = 0, 0 <= t <= 1, h = 0.2
        # actual solution: y = (1/5) t e^(3t) - (1/25) e^(3t) + (1/25) e^(-2t)
        ans1 = answer(BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2*y), 
                      0, 1, 0, 0.2,
                        Polynomial(0, 1/5) * Exponent(Polynomial(0, 3)) - (1/25) * Exponent(Polynomial(0, 3)) + (1/25) * Exponent(Polynomial(0, -2)))
        
        # (b) y' = 1 + (t - y)^2, y(2) = 1, 2 <= t <= 3, h = 0.2
        # actual solution: y = t + 1/(1-t)
        ans2 = answer(BivariateFunction(lambda t, y: 1 + (t - y)**2), 
                      2, 3, 1, 0.2,
                        Polynomial(0, 1) + 1 / (Polynomial(1, -1)))
        
        # (c) y' = 1 + y/t, y(1) = 2, 1 <= t <= 2, h = 0.2
        # actual solution: y = tlnt + 2t
        ans3 = answer(BivariateFunction(lambda t, y: 1 + y/t), 
                      1, 2, 2, 0.2,
                        Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2))
        
        return {
            'a': ans1,
            'b': ans2,
            'c': ans3
        }
    
    @staticmethod
    def problem2():
        """
        Redo problem 1 by the two step Adams Moulton implicit method.
        Compare the results with Adam-Bashforth explicit method
        """
        def answer(f, a, b, y0, h, gt):
            sol = FirstOrderLinearODE(f, a, b, y0).solve(h, method='adam-moulton', step=2, points=[gt(a+h)])(b)
            return {
                'solution': sol,
                'error': abs(gt(b) - sol)
            }
        
        # (a) y' = t e^(3t) - 2y, y(0) = 0, 0 <= t <= 1, h = 0.2
        # actual solution: y = (1/5) t e^(3t) - (1/25) e^(3t) + (1/25) e^(-2t)
        ans1 = answer(BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2*y), 
                      0, 1, 0, 0.2,
                        Polynomial(0, 1/5) * Exponent(Polynomial(0, 3)) - (1/25) * Exponent(Polynomial(0, 3)) + (1/25) * Exponent(Polynomial(0, -2)))
        
        # (b) y' = 1 + (t - y)^2, y(2) = 1, 2 <= t <= 3, h = 0.2
        # actual solution: y = t + 1/(1-t)
        ans2 = answer(BivariateFunction(lambda t, y: 1 + (t - y)**2), 
                      2, 3, 1, 0.2,
                        Polynomial(0, 1) + 1 / (Polynomial(1, -1)))
        
        # (c) y' = 1 + y/t, y(1) = 2, 1 <= t <= 2, h = 0.2
        # actual solution: y = tlnt + 2t
        ans3 = answer(BivariateFunction(lambda t, y: 1 + y/t), 
                      1, 2, 2, 0.2,
                        Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2))
        
        return {
            'a': ans1,
            'b': ans2,
            'c': ans3
        }
        
    
    @staticmethod
    def problem3():
        """
        Redo problem 1 by the three step Adams Bashforth explicit
        method and three step Adams Moulton implicit method. Compare the results.
        """
        def answer(f, a, b, y0, h, gt):
            sol1 = FirstOrderLinearODE(f, a, b, y0).solve(h, method='adam-bashforth', step=3, points=[gt(a+h), gt(a + 2 * h)])(b)
            sol2 = FirstOrderLinearODE(f, a, b, y0).solve(h, method='adam-moulton', step=3, points=[gt(a+h), gt(a + 2 * h)])(b)
            return {
                'adam-bashforth': {
                    'solution': sol1,
                    'error': abs(gt(b) - sol1)
                },
                'adam-moulton': {
                    'solution': sol2,
                    'error': abs(gt(b) - sol2)
                }
            }
        
        # (a) y' = t e^(3t) - 2y, y(0) = 0, 0 <= t <= 1, h = 0.2
        # actual solution: y = (1/5) t e^(3t) - (1/25) e^(3t) + (1/25) e^(-2t)
        ans1 = answer(BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2*y), 
                      0, 1, 0, 0.2,
                        Polynomial(0, 1/5) * Exponent(Polynomial(0, 3)) - (1/25) * Exponent(Polynomial(0, 3)) + (1/25) * Exponent(Polynomial(0, -2)))
        
        # (b) y' = 1 + (t - y)^2, y(2) = 1, 2 <= t <= 3, h = 0.2
        # actual solution: y = t + 1/(1-t)
        ans2 = answer(BivariateFunction(lambda t, y: 1 + (t - y)**2), 
                      2, 3, 1, 0.2,
                        Polynomial(0, 1) + 1 / (Polynomial(1, -1)))
        
        # (c) y' = 1 + y/t, y(1) = 2, 1 <= t <= 2, h = 0.2
        # actual solution: y = tlnt + 2t
        ans3 = answer(BivariateFunction(lambda t, y: 1 + y/t), 
                      1, 2, 2, 0.2,
                        Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2))
        
        return {
            'a': ans1,
            'b': ans2,
            'c': ans3
        }
    
    @staticmethod
    def problem4():
        """
        Apply the Adams fourth order predictor corrector method with
        h = 0.2 and starting values from the Runge Kutta fourth order
        method to the initial value problem
        y' = y - t^2 + 1, 0 <= t <= 2, y(0) = 0.5
        """
        f = BivariateFunction(lambda t, y: y - t**2 + 1)
        a = 0
        b = 2
        y0 = 0.5
        h = 0.2
        IVP = FirstOrderLinearODE(f, a, b, y0)
        sol = IVP.solve(h, method='predictor-corrector')
        
        return sol(b)
    
    @staticmethod
    def problem5():
        """
        Use the Runge Kutta method of order two to approximate the
        solution of the following problem and compare the result to the
        actual solution.
        u1' = 3u1 + 2u2 - (2t^2 + 1)e^(2t), u1(0) = 1,
        u2' = 4u1 + u2 + (t^2 + 2t - 4)e^(2t), u2(0) = 0,
        0 <= t <= 1, h = 0.2
        
        actual solution: u1 = (1/3) e^(5t) - (1/3) e^(-t) + e^(2t),
                         u2 = (1/3) e^(5t) + (2/3) e^(-t) + t^2 e^(2t)
        """
        a = 0
        b = 1
        h = 0.2
        U0 = Vector(1, 0)
        F = Vector(
            MultiVariableFunction(lambda t, u1, u2: 3*u1 + 2*u2 - (2*(t**2) + 1)*Exponent(Polynomial(0, 2))(t)),
            MultiVariableFunction(lambda t, u1, u2: 4*u1 + u2 + (t**2 + 2*t - 4)*Exponent(Polynomial(0, 2))(t))
        )
        
        IVP = FirstOrderLinearODE(F, a, b, U0)
        sol = IVP.solve(h, method='runge-kutta', n=2)
            
        GT = [
            (1/3) * Exponent(Polynomial(0, 5)) - (1/3) * Exponent(Polynomial(0, -1)) + Exponent(Polynomial(0, 2)),
            (1/3) * Exponent(Polynomial(0, 5)) + (2/3) * Exponent(Polynomial(0, -1)) + Polynomial(0, 0, 1) * Exponent(Polynomial(0, 2))
        ]
            
        return {
            'u1': {
                'result': sol[-1][0],
                'error': abs(GT[0](b) - sol[-1][0])
            },
            'u2': {
                'result': sol[-1][1],
                'error': abs(GT[1](b) - sol[-1][1])
            }
        }
