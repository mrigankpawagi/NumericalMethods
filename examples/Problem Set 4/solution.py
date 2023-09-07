from function import Function, Polynomial, Cos, Sin, Exponent
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
        Use rectangular rule and midpoint rule to evluate the integral
        int_1^5 sqrt(1+x^2) dx
        """
        f = Polynomial(1, 0, 1) ** Polynomial(0.5)
        ans1 = f.integrate(1, 5, method='rectangular')
        ans2 = f.integrate(1, 5, method='midpoint')
        return ans1, ans2

    @staticmethod
    def problem2():
        """Redo problem 1 by using the trapezoidal rule."""
        f = Polynomial(1, 0, 1) ** Polynomial(0.5)
        return f.integrate(1, 5, method='trapezoidal')

    @staticmethod
    def problem3():
        """
        Use simpson's rule to find an approximate value of the integral
        int_4^6 1 / (3 - sqrt(x)) dx
        """
        f = 1 / (3 - Polynomial(1, 0, 1) ** Polynomial(0.5))
        return f.integrate(4, 6, method='simpson')

    @staticmethod
    def problem4():
        """
        Use simpson's rule to evaluate the following.
        """
        # (a) int_0^(pi/3) cos^2(x) dx
        f = Cos(Polynomial(0, 1)) ** Polynomial(2)
        ans1 = f.integrate(0, math.pi / 3, method='simpson')

        # (b) int_0^(pi/3) sin^2(x) dx. Use answer to (a) to deduce an approximate value of this integral.
        ans2 = (math.pi / 3) - ans1

        return ans1, ans2

    @staticmethod
    def problem5():
        """
        Evaluate the integral int_0^4 (x^2 + cos x) dx by using the midpoint formula.
        """
        f = Polynomial(1, 0, 1) + Cos(Polynomial(0, 1))
        return f.integrate(0, 4, method='midpoint')

    @staticmethod
    def problem6():
        """
        Approximate the integral of f(x) = x^3 on the interval [1, 2] by using the composite trapezoidal method
        """
        f = Polynomial(0, 0, 0, 1)
        
        # (a) with four subintervals
        ans1 = f.integrate(1, 2, method='trapezoidal', n=4)

        # (b) with eight subintervals
        ans2 = f.integrate(1, 2, method='trapezoidal', n=8)

        # (c) Compute the true error in both cases
        f.integral(Polynomial(0, 0, 0, 0, 0.25))
        true_val = f.integrate(1, 2)

        return ans1, ans2, abs(ans1 - true_val), abs(ans2 - true_val)


    @staticmethod
    def problem7():
        """
        Redo problem 6 by using composite Simpson's rule.
        """
        f = Polynomial(0, 0, 0, 1)

        # (a) with four subintervals
        ans1 = f.integrate(1, 2, method='simpson', n=4)

        # (b) with eight subintervals
        ans2 = f.integrate(1, 2, method='simpson', n=8)

        # (c) Compute the true error in both cases
        f.integral(Polynomial(0, 0, 0, 0, 0.25))
        true_val = f.integrate(1, 2)

        return ans1, ans2, abs(ans1 - true_val), abs(ans2 - true_val)

    @staticmethod
    def problem8():
        """
        Using trapezoidal rule and Simpson's rule with n = 4 to approximate the value 
        of the following integral and compute the true errors and approximation errors.
        
        int_0^2 e^x^2 dx
        """
        f = Exponent(Polynomial(0, 0, 1))
        trapz = f.integrate(0, 2, method='trapezoidal', n=4)
        simps = f.integrate(0, 2, method='simpson', n=4)

        true_val = f.integrate(0, 2, 'rectangular', n=10000)

        return trapz, simps, abs(trapz - true_val), abs(simps - true_val)