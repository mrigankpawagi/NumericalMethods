from function import Function, Polynomial, Exponent, Sin, Cos, Log
import math

class Solution:

    @staticmethod
    def solve(n):
        return getattr(Solution, f'problem{n}')()

    @staticmethod
    def solve_all(n):
        return [Solution.solve(i) for i in range(1, n + 1)]

    @staticmethod
    def problem1():
        """
        f(x) = x - 2 + log(x) has a root near x = 1.5. Use the Newton-Raphson formula to obtain a better estimate.
        """
        f = Polynomial(-2, 1) + Log(Polynomial(0, 1), base=10)
        return f.root('newton', p0=1.5)
    
    @staticmethod
    def problem2():
        """
        Use Newton's method, secant method and Regula Falsi method for finding the approximations of the two zeroes, 
        one in [-1, 0] and the other in [0, 1] to withing 10^-3 accuracy of f(x) = 230x^4 + 18x^3 + 9x^2 - 221x - 9.
        Use the end points of the intervals as initial guesses for the secant method and the midpoint for Newton's 
        method.
        """
        f = Polynomial(-9, -221, 9, 18, 230)
        return ((f.root('newton', p0=-0.5, TOLERANCE=1e-3), f.root('secant', p0=-1, p1=0, TOLERANCE=1e-3), f.root('regula_falsi', p0=-1, p1=0, TOLERANCE=1e-3)), 
                (f.root('newton', p0=0.5, TOLERANCE=1e-3), f.root('secant', p0=0, p1=1, TOLERANCE=1e-3), f.root('regula_falsi', p0=0, p1=1, TOLERANCE=1e-3)))
        # The positive root is found only by Regula Falsi!

    @staticmethod
    def problem3():
        """
        Use Newton's method to find solutions accurate to within 10^-5 to the following problems.
        """
        # (a) x^3 - 2x^2 - 5 = 0 on the interval [1, 4]
        f = Polynomial(-5, 0, -2, 1)
        res1 = f.root('newton', p0=3, TOLERANCE=1e-5, return_iterations=True)

        # (b) x^2 - 2x e^-x + e^-2x = 0 on the interval [0, 1]
        f = (Polynomial(0, 1) - Exponent(Polynomial(0, -1))) ** Polynomial(2)
        # We will find the root of the derivative instead since the root of f is also the minimum and Newton's method fails around extrema.
        g = 2 * (Polynomial(0, 1) - Exponent(Polynomial(0, -1))) * (1 + Exponent(Polynomial(0, -1)))
        res2 = g.root('newton', p0=0.6, TOLERANCE=1e-5, return_iterations=True)

        # (c) x^3 - 3x^2 (2^-x) + 3x(4^-x) - 8^-x = 0 on the interval [0, 1]
        f = (Polynomial(0, 1) - Exponent(Polynomial(0, -1), base=2)) ** Polynomial(3)
        res3 = f.root('newton', p0=0.5, TOLERANCE=1e-5, return_iterations=True)

        return res1, res2, res3 

    @staticmethod
    def problem4():
        """
        Repeat the above problem using the modified Newton's method with 
        g(x) = x - f(x) * f'(x) / (f'(x)^2 - f(x) * f''(x))
        """
        # (a) x^3 - 2x^2 - 5 = 0 on the interval [1, 4]
        f = Polynomial(-5, 0, -2, 1)
        res1 = f.root('modified_newton', p0=3, TOLERANCE=1e-5, return_iterations=True)

        # (b) x^2 - 2x e^-x + e^-2x = 0 on the interval [0, 1]
        f = (Polynomial(0, 1) - Exponent(Polynomial(0, -1))) ** Polynomial(2)
        res2 = f.root('modified_newton', p0=0.6, TOLERANCE=1e-5, return_iterations=True)

        # (c) x^3 - 3x^2 (2^-x) + 3x(4^-x) - 8^-x = 0 on the interval [0, 1]
        f = (Polynomial(0, 1) - Exponent(Polynomial(0, -1), base=2)) ** Polynomial(3)
        res3 = f.root('modified_newton', p0=0.5, TOLERANCE=1e-5, return_iterations=True)

        return res1, res2, res3 
    
        # This works around minima too!
        # This takes lesser iterations on average
    
    @staticmethod
    def problem5():
        """
        Use appropriate Lagrange interpolating polynomials of degree one, two and three to approximate each of the following.
        """
        # (a) f(8.4) if f(8.1) = 16.94410, f(8.3) = 17.56492, f(8.6) = 18.50515, f(8.7) = 18.82091
        f = Polynomial.interpolate([(8.1, 16.94410), (8.3, 17.56492), (8.6, 18.50515), (8.7, 18.82091)], method='lagrange')
        res1 = f(8.4)

        # (b) f(0.25) if f(0.1) = 0.62049958, f(0.2) = -0.28398668, f(0.3) = 0.00660095, f(0.4) = 0.24842440
        f = Polynomial.interpolate([(0.1, 0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)], method='lagrange')
        res2 = f(0.25)

        return res1, res2
        # We can drop some of the points to get lower degree polynomials (these are both cubic)

    @staticmethod
    def problem6():
        """
        Construct the Lagrange interpolating polynomials for the following functions and find a bound for the absolute error
        on the interval [x_0, x_n]
        """
        # (a) f(x) = e^(2x) cos(3x), x0 = 0, x1 = 0.3, x2 = 0.6, n = 2
        f = Exponent(Polynomial(0, 2)) * Cos(Polynomial(0, 3))
        g = Polynomial.interpolate([0, 0.3, 0.6], method='lagrange', f=f)


        # (b) f(x) = sin(ln x), x0 = 2.0, x1 = 2.4, x2 = 2.6, n = 2
        f = Sin(Log(Polynomial(0, 1)))
        g = Polynomial.interpolate([2.0, 2.4, 2.6], method='lagrange', f=f)

        # (c) f(x) = ln x, x0 = 1, x1 = 1.1, x2 = 1.3, x3 = 1.4, n = 3
        f = Log(Polynomial(0, 1))
        g = Polynomial.interpolate([1, 1.1, 1.3, 1.4], method='lagrange', f=f)

        return
