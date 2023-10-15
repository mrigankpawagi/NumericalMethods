from function import Function, Polynomial, Cos, Sin, Exponent, Log, FirstOrderLinearODE
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
        # (a) y' = y/x - (y/x)^2, y(1) = 1, 1 <= x <= 2, h = 0.1
        f = lambda w: w / Polynomial(0, 1) - (w**2) / Polynomial(0, 0, 1)
        fy = lambda x: (1/x) * Polynomial(0, 1) - (1/(x**2)) * Polynomial(0, 0, 1)
        ivp = FirstOrderLinearODE(f, fy, a=1, b=2, y0=1)
        sol = ivp.solve(method='taylor', n=2, h=0.1)
        ans1 = sol(2)

    @staticmethod
    def problem2():
        pass

    @staticmethod
    def problem3():
        pass
    
    @staticmethod
    def problem4():
        pass      

    @staticmethod
    def problem5():
        pass
