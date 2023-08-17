from function import Function, Polynomial, Exponent, Sin, Cos, Tan

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
        Use the bisection method to find the solutions accurate to 1e-4 for 
        x^3 - 7x^2 + 14x - 6 = 0 on [0,1].
        """
        f = Polynomial(-6, 14, -7, 1)
        return f.root('bisection', a=0, b=1, TOLERANCE=1e-4)
    
    @staticmethod
    def problem2():
        """
        Use the bisection method to find the root of x = exp(-x) with an 
        accuracy of 1e-4. How many iterations are required?
        """
        f = Polynomial(0, 1) - Exponent(Polynomial(0, -1))
        return f.root('bisection', a=0, b=1, TOLERANCE=1e-4, return_iterations=True)

    @staticmethod
    def problem3():
        """
        Use fixed point iteration method to determine a solution accurate
        to within 1e-2 for x^3 - 3x^2 - 3 = 0 on [1,2].
        """
        f = Polynomial(-3, 0, -3, 1)
        g = Polynomial(0, 1) - f
        return g.fixed_point(p0=1.5, TOLERANCE=1e-2)

    @staticmethod
    def problem4():
        """
        Determine an interval [a,b] on which fixed point iteration will converge.
        Estimate the number of iterations required to obtain approximations 
        accurate to within 1e-5 and perform the claculations.
        """
        # (a) x = (2 - exp(x) + x^2) / 3
        # a = 0, b = 1
        f = (Polynomial(2, 0, 1) - Exponent(Polynomial(0, 1))) / Polynomial(3)
        ans_a = (0, 1), f.fixed_point(p0=0.5, TOLERANCE=1e-5)

        # (b) x = (5/x^2) + 2
        # a = 2, b = 3
        f = (Polynomial(5) / Polynomial(0, 0, 1)) + Polynomial(2)
        ans_b = (2, 3), f.fixed_point(p0=2.5, TOLERANCE=1e-5)

        # (c) x = 5^(-x)
        # a = 0, b = 1
        f = Exponent(Polynomial(0, -1), base=5)
        ans_c = (0, 1), f.fixed_point(p0=0.5, TOLERANCE=1e-5)

        # (d) x = 0.5(sinx + cosx)
        # a = 0, b = 3
        f = Polynomial(0.5) * (Sin(Polynomial(0, 1)) + Cos(Polynomial(0, 1)))
        ans_d = (0, 3), f.fixed_point(p0=1.5, TOLERANCE=1e-5)

        return ans_a, ans_b, ans_c, ans_d
    
    @staticmethod
    def problem5():
        """
        Write down the code for computing a root of a given function f(x) = 0
        using Newton Raphson's method.
        """

        # Already implemented in Function.root with method='newton'
        pass

    @staticmethod
    def problem6():
        """
        Let f(x) = -x^3 - cosx and p0 = -1. Use Newton's method to find p2.
        Could p0 = 0 be used?
        """
        f = Polynomial(0, 0, 0, -1) - Cos(Polynomial(0, 1))
        f.differentiate(Polynomial(0, 0, -3) + Sin(Polynomial(0, 1)))
        ans = f.root('newton', p0=-1, early_stop=2)

        # No, p0 = 0 cannot be used because f'(0) = 0
        return ans, False
    
    @staticmethod
    def problem7():
        """
        Use Newton's method to approximate to within 1e-4, the value of x that
        produces the point on the graph y = x^2 that is closed to (1,0).
        """
        # Distance function (squared): (x-1)^2 + (x^2-0)^2
        # We need minima of this function, that is, root of its derivative
        # which is -2 + 2x + 4x^3
        f = Polynomial(-2, 2, 0, 4)
        return f.root('newton', p0=0.5, TOLERANCE=1e-4)
    
    @staticmethod
    def problem8():
        """
        Apply Newton's method to find the approximation of the root x = tanx,
        starting with the initial guess x0 = 4 and x0 = 4.6. Compare the results
        obtained from these two intitial guesses. Does the method converge?
        """
        f = Polynomial(0, 1) - Tan(Polynomial(0, 1))
        ans_1 = f.root('newton', p0=4)
        ans_2 = f.root('newton', p0=4.6)

        return ans_1, ans_2
    
    @staticmethod
    def problem9():
        """
        Obtain an estimation (accurate till 4 decimal points) of the point of
        intersection of the curves y = x^2 - 2 and y = cosx.
        """
        f = Polynomial(-2, 0, 1) - Cos(Polynomial(0, 1))
        return f.root('newton', p0=1.5, TOLERANCE=1e-4)
    
    @staticmethod
    def problem10():
        """
        Apply Newton's method to the function f(x) = x^(2/3) if x >= 0 and
        f(x) = -x^(2/3) if x < 0, with the root x* = 0. What is the behavior of the
        iterates? Do they converge? If yes, at what order?
        """
        f = Function(lambda x: (x**(2))**(1/3) * (-1 if x < 0 else 1))
        return f.root('newton', p0=0.1, TOLERANCE=1e-5)
        # The iterates do not converge. They oscillate between two values.
