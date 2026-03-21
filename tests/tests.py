import unittest
import math

from numericalmethods import (
    Function, Polynomial, Exponent, Sin, Cos, Tan, Log,
    BivariateFunction, MultiVariableFunction,
    FirstOrderLinearODE,
    SecondOrderLinearODE_BVP, SecondOrderODE_BVP,
    Vector, Matrix, LinearSystem,
    Util,
)


class TestProblemSet1(unittest.TestCase):
    """Tests for Problem Set 1 - Root Finding Methods"""

    def test_problem1(self):
        """
        Use the bisection method to find the solutions accurate to 1e-4 for
        x^3 - 7x^2 + 14x - 6 = 0 on [0, 1].
        """
        f = Polynomial(-6, 14, -7, 1)
        result = f.root('bisection', a=0, b=1, TOLERANCE=1e-4)
        self.assertAlmostEqual(result, 0.585784912109375, places=10)

    def test_problem2(self):
        """
        Use the bisection method to find the root of x = exp(-x) with an
        accuracy of 1e-4. How many iterations are required?
        """
        f = Polynomial(0, 1) - Exponent(Polynomial(0, -1))
        root, iterations = f.root('bisection', a=0, b=1, TOLERANCE=1e-4, return_iterations=True)
        self.assertAlmostEqual(root, 0.567169189453125, places=10)
        self.assertEqual(iterations, 15)

    def test_problem3(self):
        """
        Use fixed point iteration method to determine a solution accurate
        to within 1e-2 for x^3 - 3x^2 - 3 = 0 on [1, 2].
        The iteration does not converge for this setup.
        """
        f = Polynomial(-3, 0, -3, 1)
        g = Polynomial(0, 1) - f
        result = g.fixed_point(p0=1.5, TOLERANCE=1e-2)
        self.assertIsNone(result)

    def test_problem4(self):
        """
        Determine an interval [a, b] on which fixed point iteration will converge.
        Estimate the number of iterations required to obtain approximations
        accurate to within 1e-5 and perform the calculations.
        (a) x = (2 - exp(x) + x^2) / 3, on [0, 1]
        (b) x = (5 / x^2) + 2, on [2, 3]
        (c) x = 5^(-x), on [0, 1]
        (d) x = 0.5 * (sin x + cos x), on [0, 3]
        """
        # (a) x = (2 - exp(x) + x^2) / 3
        f = (Polynomial(2, 0, 1) - Exponent(Polynomial(0, 1))) / Polynomial(3)
        result_a = f.fixed_point(p0=0.5, TOLERANCE=1e-5)
        self.assertAlmostEqual(result_a, 0.25752908416795584, delta=1e-5)

        # (b) x = 5 / x^2 + 2
        f = (Polynomial(5) / Polynomial(0, 0, 1)) + Polynomial(2)
        result_b = f.fixed_point(p0=2.5, TOLERANCE=1e-5)
        self.assertAlmostEqual(result_b, 2.6906498932822, delta=1e-5)

        # (c) x = 5^(-x)
        f = Exponent(Polynomial(0, -1), base=5)
        result_c = f.fixed_point(p0=0.5, TOLERANCE=1e-5)
        self.assertAlmostEqual(result_c, 0.4696257766724794, delta=1e-5)

        # (d) x = 0.5 * (sin x + cos x)
        f = Polynomial(0.5) * (Sin(Polynomial(0, 1)) + Cos(Polynomial(0, 1)))
        result_d = f.fixed_point(p0=1.5, TOLERANCE=1e-5)
        self.assertAlmostEqual(result_d, 0.7048117652456208, delta=1e-5)

    def test_problem6(self):
        """
        Let f(x) = -x^3 - cos x and p0 = -1. Use Newton's method to find p2.
        Could p0 = 0 be used? (No, because f'(0) = 0.)
        """
        f = Polynomial(0, 0, 0, -1) - Cos(Polynomial(0, 1))
        f.differentiate(Polynomial(0, 0, -3) + Sin(Polynomial(0, 1)))
        root = f.root('newton', p0=-1, early_stop=2)
        self.assertAlmostEqual(root, -0.8654740738468996, places=10)

    def test_problem7(self):
        """
        Use Newton's method to approximate to within 1e-4 the value of x that
        produces the point on the graph y = x^2 closest to (1, 0).
        Minimize the squared distance (x-1)^2 + x^4, i.e. find the root of
        its derivative -2 + 2x + 4x^3.
        """
        f = Polynomial(-2, 2, 0, 4)
        result = f.root('newton', p0=0.5, TOLERANCE=1e-4)
        self.assertAlmostEqual(result, 0.5897545123016604, delta=1e-4)

    def test_problem8(self):
        """
        Apply Newton's method to find the approximation of the root x = tan x,
        starting with x0 = 4 and x0 = 4.6.
        Starting from x0 = 4 Newton's method does not converge; x0 = 4.6 converges.
        """
        f = Polynomial(0, 1) - Tan(Polynomial(0, 1))
        result_1 = f.root('newton', p0=4)
        result_2 = f.root('newton', p0=4.6)
        self.assertIsNone(result_1)
        self.assertAlmostEqual(result_2, 4.493409457909064, places=7)

    def test_problem9(self):
        """
        Obtain an estimation (accurate to 4 decimal points) of the point of
        intersection of the curves y = x^2 - 2 and y = cos x.
        """
        f = Polynomial(-2, 0, 1) - Cos(Polynomial(0, 1))
        result = f.root('newton', p0=1.5, TOLERANCE=1e-4)
        self.assertAlmostEqual(result, 1.4546189292083291, delta=1e-4)

    def test_problem10(self):
        """
        Apply Newton's method to f(x) = x^(2/3) * sign(x) with root x* = 0.
        The method converges to near zero (at order 1).
        """
        f = Function(lambda x: (x ** 2) ** (1 / 3) * (-1 if x < 0 else 1))
        result = f.root('newton', p0=0.1, TOLERANCE=1e-5)
        self.assertAlmostEqual(result, 0.0, delta=1e-5)


class TestProblemSet2(unittest.TestCase):
    """Tests for Problem Set 2 - Newton's Method, Secant, Regula Falsi, Interpolation"""

    def test_problem1(self):
        """
        f(x) = x - 2 + log10(x) has a root near x = 1.5.
        Use the Newton-Raphson formula to obtain a better estimate.
        """
        f = Polynomial(-2, 1) + Log(Polynomial(0, 1), base=10)
        result = f.root('newton', p0=1.5)
        self.assertAlmostEqual(result, 1.7555794992611777, places=7)

    def test_problem2(self):
        """
        Use Newton's method, secant method, and Regula Falsi to find the two zeroes
        of f(x) = 230x^4 + 18x^3 + 9x^2 - 221x - 9 in [-1, 0] and [0, 1].
        Use tolerance 1e-3.
        """
        f = Polynomial(-9, -221, 9, 18, 230)

        # Root in [-1, 0]
        newton_neg = f.root('newton', p0=-0.5, TOLERANCE=1e-3)
        secant_neg = f.root('secant', p0=-1, p1=0, TOLERANCE=1e-3)
        rf_neg = f.root('regula_falsi', p0=-1, p1=0, TOLERANCE=1e-3)
        self.assertAlmostEqual(newton_neg, -0.04065928831573657, places=3)
        self.assertAlmostEqual(secant_neg, -0.04065926257769109, places=3)
        self.assertAlmostEqual(rf_neg, -0.0399800081852857, places=3)

        # Root in [0, 1] — only Regula Falsi converges to the positive root
        newton_pos = f.root('newton', p0=0.5, TOLERANCE=1e-3)
        secant_pos = f.root('secant', p0=0, p1=1, TOLERANCE=1e-3)
        rf_pos = f.root('regula_falsi', p0=0, p1=1, TOLERANCE=1e-3)
        self.assertAlmostEqual(newton_pos, -0.04065928833435429, places=3)
        self.assertAlmostEqual(secant_pos, -0.04065922824320606, places=3)
        self.assertAlmostEqual(rf_pos, 0.962391746837012, places=3)

    def test_problem3(self):
        """
        Use Newton's method to find solutions accurate to within 1e-5.
        (a) x^3 - 2x^2 - 5 = 0 on [1, 4]
        (b) (x - e^{-x})^2 = 0 on [0, 1]  (use derivative to avoid the extremum)
        (c) (x - 2^{-x})^3 = 0 on [0, 1]
        """
        # (a)
        f = Polynomial(-5, 0, -2, 1)
        root_a, iters_a = f.root('newton', p0=3, TOLERANCE=1e-5, return_iterations=True)
        self.assertAlmostEqual(root_a, 2.690647448031735, places=5)
        self.assertEqual(iters_a, 4)

        # (b) solve via the derivative of f to avoid issues around the minimum
        g = 2 * (Polynomial(0, 1) - Exponent(Polynomial(0, -1))) * (1 + Exponent(Polynomial(0, -1)))
        root_b, iters_b = g.root('newton', p0=0.6, TOLERANCE=1e-5, return_iterations=True)
        self.assertAlmostEqual(root_b, 0.5671432904107577, places=5)
        self.assertEqual(iters_b, 3)

        # (c)
        f = (Polynomial(0, 1) - Exponent(Polynomial(0, -1), base=2)) ** Polynomial(3)
        root_c, iters_c = f.root('newton', p0=0.5, TOLERANCE=1e-5, return_iterations=True)
        self.assertAlmostEqual(root_c, 0.641187997338787, places=5)
        self.assertEqual(iters_c, 23)

    def test_problem4(self):
        """
        Repeat Problem 3 using the modified Newton's method
        g(x) = x - f(x)*f'(x) / (f'(x)^2 - f(x)*f''(x)).
        This method also works near minima and typically takes fewer iterations.
        """
        # (a)
        f = Polynomial(-5, 0, -2, 1)
        root_a, iters_a = f.root('modified_newton', p0=3, TOLERANCE=1e-5, return_iterations=True)
        self.assertAlmostEqual(root_a, 2.6906474480239617, places=5)
        self.assertEqual(iters_a, 4)

        # (b)
        f = (Polynomial(0, 1) - Exponent(Polynomial(0, -1))) ** Polynomial(2)
        root_b, iters_b = f.root('modified_newton', p0=0.6, TOLERANCE=1e-5, return_iterations=True)
        self.assertAlmostEqual(root_b, 0.5671501775692881, places=5)
        self.assertEqual(iters_b, 3)

        # (c)
        f = (Polynomial(0, 1) - Exponent(Polynomial(0, -1), base=2)) ** Polynomial(3)
        root_c, iters_c = f.root('modified_newton', p0=0.5, TOLERANCE=1e-5, return_iterations=True)
        self.assertAlmostEqual(root_c, 0.641197156792501, places=5)
        self.assertEqual(iters_c, 4)

    def test_problem5(self):
        """
        Use appropriate Lagrange interpolating polynomials of degree 1, 2, and 3
        to approximate the following.
        (a) f(8.4) from data at x = 8.1, 8.3, 8.6, 8.7
        (b) f(0.25) from data at x = 0.1, 0.2, 0.3, 0.4
        """
        # (a) cubic Lagrange through four points
        f = Polynomial.interpolate(
            [(8.1, 16.94410), (8.3, 17.56492), (8.6, 18.50515), (8.7, 18.82091)],
            method='lagrange'
        )
        self.assertAlmostEqual(f(8.4), 17.8771425, places=6)

        # (b) cubic Lagrange through four points
        f = Polynomial.interpolate(
            [(0.1, 0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)],
            method='lagrange'
        )
        self.assertAlmostEqual(f(0.25), -0.21033722187499995, places=6)


class TestProblemSet3(unittest.TestCase):
    """Tests for Problem Set 3 - Forward/Backward Difference and Lagrange Interpolation"""

    def test_problem1(self):
        """
        Use Newton's forward difference formula to construct interpolating polynomials
        of degree 1, 2, and 3 for the following data.
        (a) Estimate f(-1/3) from data at x = -0.75, -0.5, -0.25, 0.
        (b) Estimate f(0.25) from data at x = 0.1, 0.2, 0.3, 0.4.
        """
        # (a)
        data_a = [(-0.75, -0.07181250), (-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101000)]
        p1 = Polynomial.interpolate(data_a[:2], form='forward_diff')
        p2 = Polynomial.interpolate(data_a[:3], form='forward_diff')
        p3 = Polynomial.interpolate(data_a[:4], form='forward_diff')
        self.assertAlmostEqual(p1(-1 / 3), 0.006625000000000006, places=10)
        self.assertAlmostEqual(p2(-1 / 3), 0.1803055555555556, places=10)
        self.assertAlmostEqual(p3(-1 / 3), 0.17451851851851857, places=10)

        # (b)
        data_b = [(0.1, -0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)]
        p1 = Polynomial.interpolate(data_b[:2], form='forward_diff')
        p2 = Polynomial.interpolate(data_b[:3], form='forward_diff')
        p3 = Polynomial.interpolate(data_b[:4], form='forward_diff')
        self.assertAlmostEqual(p1(0.25), -0.11573022999999993, places=10)
        self.assertAlmostEqual(p2(0.25), -0.13295220624999998, places=10)
        self.assertAlmostEqual(p3(0.25), -0.132774774375, places=10)

    def test_problem2(self):
        """
        Redo Problem 1 using Newton's backward difference formula.
        """
        # (a)
        data_a = [(-0.75, -0.07181250), (-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101000)]
        p1 = Polynomial.interpolate(data_a[:2], form='backward_diff')
        p2 = Polynomial.interpolate(data_a[:3], form='backward_diff')
        p3 = Polynomial.interpolate(data_a[:4], form='backward_diff')
        self.assertAlmostEqual(p1(-1 / 3), 0.006625000000000006, places=10)
        self.assertAlmostEqual(p2(-1 / 3), 0.18030555555555558, places=10)
        self.assertAlmostEqual(p3(-1 / 3), 0.1745185185185185, places=10)

        # (b)
        data_b = [(0.1, -0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)]
        p1 = Polynomial.interpolate(data_b[:2], form='backward_diff')
        p2 = Polynomial.interpolate(data_b[:3], form='backward_diff')
        p3 = Polynomial.interpolate(data_b[:4], form='backward_diff')
        self.assertAlmostEqual(p1(0.25), -0.11573022999999996, places=10)
        self.assertAlmostEqual(p2(0.25), -0.13295220624999998, places=10)
        self.assertAlmostEqual(p3(0.25), -0.1327747743750001, places=10)

    def test_problem3(self):
        """
        Find the degree of the polynomial that interpolates
        f(-2)=1, f(-1)=4, f(0)=11, f(1)=16, f(2)=13, f(3)=-4.
        The forward differences reveal that the polynomial is of degree 3.
        """
        data = [(-2, 1), (-1, 4), (0, 11), (1, 16), (2, 13), (3, -4)]
        self.assertEqual(Util.delta(1, 0, data), 3)
        self.assertEqual(Util.delta(2, 0, data), 4)
        self.assertEqual(Util.delta(3, 0, data), -6)
        self.assertEqual(Util.delta(4, 0, data), 0)
        self.assertEqual(Util.delta(5, 0, data), 0)

    def test_problem4(self):
        """
        Use Lagrange interpolating polynomials of degree 1, 2, and 3 to
        approximate the data given in Problem 1.
        (a) Estimate f(-1/3) from data at x = -0.75, -0.5, -0.25, 0.
        (b) Estimate f(0.25) from data at x = 0.1, 0.2, 0.3, 0.4.
        """
        # (a)
        data_a = [(-0.75, -0.07181250), (-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101000)]
        p1 = Polynomial.interpolate(data_a[:2], method='lagrange')
        p2 = Polynomial.interpolate(data_a[:3], method='lagrange')
        p3 = Polynomial.interpolate(data_a[:4], method='lagrange')
        self.assertAlmostEqual(p1(-1 / 3), 0.006625000000000006, places=10)
        self.assertAlmostEqual(p2(-1 / 3), 0.1803055555555556, places=10)
        self.assertAlmostEqual(p3(-1 / 3), 0.1745185185185186, places=10)

        # (b)
        data_b = [(0.1, -0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)]
        p1 = Polynomial.interpolate(data_b[:2], method='lagrange')
        p2 = Polynomial.interpolate(data_b[:3], method='lagrange')
        p3 = Polynomial.interpolate(data_b[:4], method='lagrange')
        self.assertAlmostEqual(p1(0.25), -0.11573022999999993, places=10)
        self.assertAlmostEqual(p2(0.25), -0.13295220624999998, places=10)
        self.assertAlmostEqual(p3(0.25), -0.132774774375, places=10)


class TestProblemSet4(unittest.TestCase):
    """Tests for Problem Set 4 - Numerical Integration"""

    def test_problem1(self):
        """
        Use rectangular rule and midpoint rule to evaluate
        int_1^5 sqrt(1 + x^2) dx.
        """
        f = Polynomial(1, 0, 1) ** Polynomial(0.5)
        rect = f.integrate(1, 5, method='rectangular')
        mid = f.integrate(1, 5, method='midpoint')
        self.assertAlmostEqual(rect, 5.656854249492381, places=7)
        self.assertAlmostEqual(mid, 12.649110640673518, places=7)

    def test_problem2(self):
        """
        Redo Problem 1 using the trapezoidal rule.
        """
        f = Polynomial(1, 0, 1) ** Polynomial(0.5)
        result = f.integrate(1, 5, method='trapezoidal')
        self.assertAlmostEqual(result, 13.026466151931759, places=7)

    def test_problem3(self):
        """
        Use Simpson's rule to find an approximate value of
        int_4^6 1 / (3 - sqrt(x)) dx.
        """
        f = 1 / (3 - Polynomial(0, 1) ** Polynomial(0.5))
        result = f.integrate(4, 6, method='simpson')
        self.assertAlmostEqual(result, 2.684188186142505, places=7)

    def test_problem4(self):
        """
        Use Simpson's rule to evaluate:
        (a) int_0^{pi/3} cos^2(x) dx
        (b) int_0^{pi/3} sin^2(x) dx  (deduced using cos^2 + sin^2 = 1)
        """
        f = Cos(Polynomial(0, 1)) ** Polynomial(2)
        ans1 = f.integrate(0, math.pi / 3, method='simpson')
        ans2 = (math.pi / 3) - ans1
        self.assertAlmostEqual(ans1, 0.74176493209759, places=7)
        self.assertAlmostEqual(ans2, 0.30543261909900765, places=7)

    def test_problem5(self):
        """
        Evaluate int_0^4 (x^2 + cos x) dx using the midpoint formula.
        """
        f = Polynomial(0, 0, 1) + Cos(Polynomial(0, 1))
        result = f.integrate(0, 4, method='midpoint')
        self.assertAlmostEqual(result, 14.33541265381143, places=7)

    def test_problem6(self):
        """
        Approximate int_1^2 x^3 dx using the composite trapezoidal method.
        (a) with four subintervals
        (b) with eight subintervals
        (c) Compute the true error in both cases.  True value = 3.75.
        """
        f = Polynomial(0, 0, 0, 1)
        ans1 = f.integrate(1, 2, method='trapezoidal', n=4)
        ans2 = f.integrate(1, 2, method='trapezoidal', n=8)
        f.integral(Polynomial(0, 0, 0, 0, 0.25))
        true_val = f.integrate(1, 2)
        self.assertAlmostEqual(ans1, 3.796875, places=7)
        self.assertAlmostEqual(ans2, 3.76171875, places=7)
        self.assertAlmostEqual(abs(ans1 - true_val), 0.046875, places=7)
        self.assertAlmostEqual(abs(ans2 - true_val), 0.01171875, places=7)

    def test_problem7(self):
        """
        Redo Problem 6 using composite Simpson's rule.
        (a) with four subintervals
        (b) with eight subintervals
        (c) Compute the true error in both cases.  Simpson's rule is exact for cubics.
        """
        f = Polynomial(0, 0, 0, 1)
        ans1 = f.integrate(1, 2, method='simpson', n=4)
        ans2 = f.integrate(1, 2, method='simpson', n=8)
        f.integral(Polynomial(0, 0, 0, 0, 0.25))
        true_val = f.integrate(1, 2)
        self.assertAlmostEqual(ans1, 3.75, places=7)
        self.assertAlmostEqual(ans2, 3.75, places=7)
        self.assertAlmostEqual(abs(ans1 - true_val), 0.0, places=7)
        self.assertAlmostEqual(abs(ans2 - true_val), 0.0, places=7)

    def test_problem8(self):
        """
        Use trapezoidal rule and Simpson's rule with n = 4 to approximate
        int_0^2 e^{x^2} dx, and compute the true errors.
        """
        f = Exponent(Polynomial(0, 0, 1))
        trapz = f.integrate(0, 2, method='trapezoidal', n=4)
        simps = f.integrate(0, 2, method='simpson', n=4)
        true_val = f.integrate(0, 2, 'rectangular', n=10000)
        self.assertAlmostEqual(trapz, 20.644559049038712, places=7)
        self.assertAlmostEqual(simps, 16.538594702002605, places=7)
        self.assertAlmostEqual(abs(trapz - true_val), 4.197290370559472, places=3)
        self.assertAlmostEqual(abs(simps - true_val), 0.09132602352336505, places=3)


class TestProblemSet5(unittest.TestCase):
    """Tests for Problem Set 5 - Composite Integration and Gaussian Quadrature"""

    def test_problem1(self):
        """
        Approximate int_1^5 (x^3 + 5x^2 + 1) dx using composite rectangular
        and composite midpoint methods with n = 5 and n = 10 subintervals.
        Also compute the true errors.
        """
        f = Polynomial(1, 0, 5, 1)
        res_rect_5 = f.integrate(1, 5, method='rectangular', n=5)
        res_mid_5 = f.integrate(1, 5, method='midpoint', n=5)
        res_rect_10 = f.integrate(1, 5, method='rectangular', n=10)
        res_mid_10 = f.integrate(1, 5, method='midpoint', n=10)
        f.integral(Polynomial(0, 1, 0, 5 / 3, 1 / 4))
        true = f.integrate(1, 5)

        self.assertAlmostEqual(res_rect_5, 275.04, places=7)
        self.assertAlmostEqual(res_mid_5, 363.68000000000006, places=7)
        self.assertAlmostEqual(res_rect_10, 319.36, places=7)
        self.assertAlmostEqual(res_mid_10, 365.9200000000001, places=7)
        self.assertAlmostEqual(abs(true - res_rect_5), 91.62666666666667, places=5)
        self.assertAlmostEqual(abs(true - res_mid_5), 2.986666666666622, places=5)
        self.assertAlmostEqual(abs(true - res_rect_10), 47.30666666666667, places=5)
        self.assertAlmostEqual(abs(true - res_mid_10), 0.7466666666666129, places=5)

    def test_problem2(self):
        """
        Redo Problem 1 using composite trapezoidal and composite Simpson methods.
        Simpson's rule integrates cubics exactly, so the true error is 0.
        """
        f = Polynomial(1, 0, 5, 1)
        res_trap_5 = f.integrate(1, 5, method='trapezoidal', n=5)
        res_simp_5 = f.integrate(1, 5, method='simpson', n=5)
        res_trap_10 = f.integrate(1, 5, method='trapezoidal', n=10)
        res_simp_10 = f.integrate(1, 5, method='simpson', n=10)
        f.integral(Polynomial(0, 1, 0, 5 / 3, 1 / 4))
        true = f.integrate(1, 5)

        self.assertAlmostEqual(res_trap_5, 372.64000000000004, places=7)
        self.assertAlmostEqual(res_simp_5, 366.6666666666667, places=7)
        self.assertAlmostEqual(res_trap_10, 368.16, places=7)
        self.assertAlmostEqual(res_simp_10, 366.6666666666667, places=7)
        self.assertAlmostEqual(abs(true - res_trap_5), 5.973333333333358, places=5)
        self.assertAlmostEqual(abs(true - res_simp_5), 0.0, places=7)
        self.assertAlmostEqual(abs(true - res_trap_10), 1.4933333333333394, places=5)
        self.assertAlmostEqual(abs(true - res_simp_10), 0.0, places=7)

    def test_problem3(self):
        """
        Evaluate int_0^{pi/2} x sin(x) dx using one-point Gauss quadrature,
        and compute the true error.  True value = 1 (via integration by parts).
        """
        f = Polynomial(0, 1) * Sin(Polynomial(0, 1))
        f.integral(Sin(Polynomial(0, 1)) - Polynomial(0, 1) * Cos(Polynomial(0, 1)))
        res = f.integrate(0, math.pi / 2, method='gauss', n=1)
        true = f.integrate(0, math.pi / 2)
        self.assertAlmostEqual(res, 0.8723580249548599, places=7)
        self.assertAlmostEqual(abs(true - res), 0.12764197504513997, places=7)

    def test_problem4(self):
        """
        Redo Problem 3 using two-point Gauss quadrature formula.
        """
        f = Polynomial(0, 1) * Sin(Polynomial(0, 1))
        f.integral(Sin(Polynomial(0, 1)) - Polynomial(0, 1) * Cos(Polynomial(0, 1)))
        res = f.integrate(0, math.pi / 2, method='gauss', n=2)
        true = f.integrate(0, math.pi / 2)
        self.assertAlmostEqual(res, 1.0048348693320484, places=7)
        self.assertAlmostEqual(abs(true - res), 0.004834869332048464, places=7)

    def test_problem5(self):
        """
        Use composite Simpson's rule with m = 2, n = 4 to approximate
        int_{1.4}^{2.0} int_{1.0}^{1.5} ln(x + 2y) dy dx.
        """
        m, n = 2, 4
        a, b, c, d = 1.4, 2, 1, 1.5
        h = (d - c) / m
        f = lambda r: Log(Polynomial(2 * r, 1))
        g = (h / 6) * (
            f(d) + f(c)
            + 2 * sum(f(c + i * h) for i in range(1, m))
            + 4 * sum(f(c + i * h + h / 2) for i in range(m))
        )
        res = g.integrate(a, b, method='simpson', n=n)
        self.assertAlmostEqual(res, 0.42955439362739267, places=7)

    def test_problem6(self):
        """
        Redo Problem 5 using the Gaussian quadrature formula with n = 1
        in both dimensions.
        """
        a, b = 1.0, 1.5
        c, d = 1.4, 2.0
        f = lambda r: ((b - a) / 2) * Log(
            Polynomial(2 * (((a + b) / 2) + ((b - a) / 2) * r), 1)
        )
        g = f(-1 / math.sqrt(3)) + f(1 / math.sqrt(3))
        res = g.integrate(c, d, method='gauss', n=1)
        self.assertAlmostEqual(res, 0.42981506172407835, places=7)


class TestProblemSet6(unittest.TestCase):
    """Tests for Problem Set 6 - Numerical Differentiation and Euler's Method"""

    def test_problem1(self):
        """
        Use the forward difference formula to approximate the derivative of
        f(x) = ln x at x0 = 1.8 using h = 0.1, 0.05, and 0.01.
        Also determine the error bounds.
        """
        f = Log(Polynomial(0, 1))
        a = 1.8
        max_err = lambda h: h / (2 * a ** 2)

        res_01 = f.differentiate(h=0.1, method='forward')(a)
        res_005 = f.differentiate(h=0.05, method='forward')(a)
        res_001 = f.differentiate(h=0.01, method='forward')(a)

        self.assertAlmostEqual(res_01, 0.5406722127027574, places=7)
        self.assertAlmostEqual(res_005, 0.5479794837622887, places=7)
        self.assertAlmostEqual(res_001, 0.5540180375615322, places=7)
        self.assertAlmostEqual(max_err(0.1), 0.015432098765432098, places=7)
        self.assertAlmostEqual(max_err(0.05), 0.007716049382716049, places=7)
        self.assertAlmostEqual(max_err(0.01), 0.0015432098765432098, places=7)

    def test_problem2(self):
        """
        Redo Problem 1 with the backward difference formula and the
        central difference formula.
        """
        f = Log(Polynomial(0, 1))
        a = 1.8

        # backward difference
        back_max_err = lambda h: h / (2 * (a - h) ** 2)
        res_b01 = f.differentiate(h=0.1, method='backward')(a)
        res_b005 = f.differentiate(h=0.05, method='backward')(a)
        res_b001 = f.differentiate(h=0.01, method='backward')(a)
        self.assertAlmostEqual(res_b01, 0.5715841383994869, places=7)
        self.assertAlmostEqual(res_b005, 0.563417539333928, places=7)
        self.assertAlmostEqual(res_b001, 0.5571045049455381, places=7)
        self.assertAlmostEqual(back_max_err(0.1), 0.01730103806228374, places=7)
        self.assertAlmostEqual(back_max_err(0.05), 0.00816326530612245, places=7)
        self.assertAlmostEqual(back_max_err(0.01), 0.0015605006085952374, places=7)

        # central difference
        cent_max_err = lambda h: (h ** 2) / (3 * (a - h) ** 3)
        res_c01 = f.differentiate(h=0.1, method='central')(a)
        res_c005 = f.differentiate(h=0.05, method='central')(a)
        res_c001 = f.differentiate(h=0.01, method='central')(a)
        self.assertAlmostEqual(res_c01, 0.5561281755511222, places=7)
        self.assertAlmostEqual(res_c005, 0.5556985115481083, places=7)
        self.assertAlmostEqual(res_c001, 0.5555612712535352, places=7)
        self.assertAlmostEqual(cent_max_err(0.1), 0.0006784720808738723, places=7)
        self.assertAlmostEqual(cent_max_err(0.05), 0.00015549076773566572, places=7)
        self.assertAlmostEqual(cent_max_err(0.01), 5.8119203299636395e-06, places=15)

    def test_problem3(self):
        """
        Consider the IVP y' = y * ln y / x, y(2) = e.
        Use Euler's method with h = 0.1 to approximate y(3).
        """
        f = BivariateFunction(
            lambda x, y: Polynomial(0, 1)(y) * Log(Polynomial(0, 1))(y) / Polynomial(0, 1)(x)
        )
        IVP = FirstOrderLinearODE(f, 2, 3, math.e)
        sol = IVP.solve(h=0.1, method='euler')
        self.assertAlmostEqual(sol(3), 4.418072257635635, places=5)

    def test_problem4(self):
        """
        Consider the IVP y' = y - x, y(0) = 1/2.
        Use Euler's method with h = 0.1 and h = 0.05 to approximate y(1).
        The exact solution is y(x) = x + 1 - (1/2) e^x.
        """
        f = BivariateFunction(lambda x, y: y - x)
        gt = (Polynomial(1, 1) - 0.5 * Exponent(Polynomial(0, 1)))(1)
        IVP = FirstOrderLinearODE(f, 0, 1, 0.5)

        sol1 = IVP.solve(h=0.1, method='euler')(1)
        sol2 = IVP.solve(h=0.05, method='euler')(1)

        self.assertAlmostEqual(sol1, 0.7031287699500001, places=7)
        self.assertAlmostEqual(abs(gt - sol1), 0.06226968417952261, places=5)
        self.assertAlmostEqual(sol2, 0.6733511474277903, places=7)
        self.assertAlmostEqual(abs(gt - sol2), 0.03249206165731289, places=5)


class TestProblemSet7(unittest.TestCase):
    """Tests for Problem Set 7 - Taylor's Method, Runge-Kutta, and Trapezoidal ODE Solvers"""

    def _solve(self, f, a, b, y0, h, method, n=None):
        ode = FirstOrderLinearODE(f, a, b, y0)
        if n is not None:
            return ode.solve(h, method=method, n=n)(b)
        return ode.solve(h, method=method)(b)

    def test_problem1(self):
        """
        Use Taylor's series method of order 2 to approximate the solution for
        each of the following IVPs.
        (a) y' = y/x - (y/x)^2, y(1) = 1, h = 0.1
        (b) y' = sin x + e^{-x}, y(0) = 0, h = 0.5
        (c) y' = (y^2 + y) / x, y(1) = -2, h = 0.5
        (d) y' = -xy + 4x/y, y(0) = 1, h = 0.25
        """
        ans1 = self._solve(
            BivariateFunction(lambda x, y: y / x - (y / x) ** 2), 1, 2, 1, 0.1, 'taylor', n=2
        )
        ans2 = self._solve(
            BivariateFunction(lambda x, y: Sin(Polynomial(0, 1))(x) + Exponent(Polynomial(0, -1))(x)),
            0, 1, 0, 0.5, 'taylor', n=2
        )
        ans3 = self._solve(
            BivariateFunction(lambda x, y: (y ** 2 + y) / x), 1, 3, -2, 0.5, 'taylor', n=2
        )
        ans4 = self._solve(
            BivariateFunction(lambda x, y: -x * y + 4 * x / y), 0, 1, 1, 0.25, 'taylor', n=2
        )
        self.assertAlmostEqual(ans1, 1.1827423857604942, places=7)
        self.assertAlmostEqual(ans2, 1.0768602913630998, places=7)
        self.assertAlmostEqual(ans3, -1.4588804494524443, places=7)
        self.assertAlmostEqual(ans4, 1.7204836614121437, places=7)

    def test_problem2(self):
        """
        Redo Problem 1 using the Runge-Kutta method of order 2.
        """
        ans1 = self._solve(
            BivariateFunction(lambda x, y: y / x - (y / x) ** 2), 1, 2, 1, 0.1, 'runge-kutta', n=2
        )
        ans2 = self._solve(
            BivariateFunction(lambda x, y: Sin(Polynomial(0, 1))(x) + Exponent(Polynomial(0, -1))(x)),
            0, 1, 0, 0.5, 'runge-kutta', n=2
        )
        ans3 = self._solve(
            BivariateFunction(lambda x, y: (y ** 2 + y) / x), 1, 3, -2, 0.5, 'runge-kutta', n=2
        )
        ans4 = self._solve(
            BivariateFunction(lambda x, y: -x * y + 4 * x / y), 0, 1, 1, 0.25, 'runge-kutta', n=2
        )
        self.assertAlmostEqual(ans1, 1.1808344690528974, places=7)
        self.assertAlmostEqual(ans2, 1.0953157056532528, places=7)
        self.assertAlmostEqual(ans3, -1.2020872542290026, places=7)
        self.assertAlmostEqual(ans4, 1.6922889668511303, places=7)

    def test_problem3(self):
        """
        Redo Problem 1 using the Trapezoidal method.
        """
        ans1 = self._solve(
            BivariateFunction(lambda x, y: y / x - (y / x) ** 2), 1, 2, 1, 0.1, 'trapezoidal'
        )
        ans2 = self._solve(
            BivariateFunction(lambda x, y: Sin(Polynomial(0, 1))(x) + Exponent(Polynomial(0, -1))(x)),
            0, 1, 0, 0.5, 'trapezoidal'
        )
        ans3 = self._solve(
            BivariateFunction(lambda x, y: (y ** 2 + y) / x), 1, 3, -2, 0.5, 'trapezoidal'
        )
        ans4 = self._solve(
            BivariateFunction(lambda x, y: -x * y + 4 * x / y), 0, 1, 1, 0.25, 'trapezoidal'
        )
        self.assertAlmostEqual(ans1, 1.1805847365032, places=7)
        self.assertAlmostEqual(ans2, 1.095315705653253, places=7)
        self.assertAlmostEqual(ans3, -1.1614093052592742, places=7)
        self.assertAlmostEqual(ans4, 1.6948357863442574, places=7)

    def test_problem4(self):
        """
        Compare Taylor's (order 2), Runge-Kutta (order 2), and Trapezoidal methods
        on three IVPs with known exact solutions.
        (a) y' = x e^{3x} - 2y, exact: y = (1/5)x e^{3x} - (1/25)e^{3x} + (1/25)e^{-2x}
        (b) y' = 1 + (x - y)^2, exact: y = x + 1/(1 - x)
        (c) y' = 1 + y/x, exact: y = x ln x + 2x
        """
        # (a)
        fa = BivariateFunction(lambda x, y: x * Exponent(Polynomial(0, 3))(x) - 2 * y)
        gt_a = (
            Polynomial(0, 1 / 5) * Exponent(Polynomial(0, 3))
            - (1 / 25) * Exponent(Polynomial(0, 3))
            + (1 / 25) * Exponent(Polynomial(0, -2))
        )
        sol_a_t = self._solve(fa, 0, 1, 0, 0.1, 'taylor', n=2)
        sol_a_r = self._solve(fa, 0, 1, 0, 0.1, 'runge-kutta', n=2)
        sol_a_p = self._solve(fa, 0, 1, 0, 0.1, 'trapezoidal')
        self.assertAlmostEqual(sol_a_t, 3.161454578035206, places=7)
        self.assertAlmostEqual(sol_a_r, 3.297890507632929, places=7)
        self.assertAlmostEqual(sol_a_p, 3.2477418088160586, places=7)
        self.assertAlmostEqual(abs(gt_a(1) - sol_a_t), 0.05764474100428485, places=5)
        self.assertAlmostEqual(abs(gt_a(1) - sol_a_r), 0.07879118859343803, places=5)
        self.assertAlmostEqual(abs(gt_a(1) - sol_a_p), 0.028642489776567803, places=5)

        # (b)
        fb = BivariateFunction(lambda x, y: 1 + (x - y) ** 2)
        gt_b = Polynomial(0, 1) + 1 / Polynomial(1, -1)
        sol_b_t = self._solve(fb, 2, 3, 1, 0.5, 'taylor', n=2)
        sol_b_r = self._solve(fb, 2, 3, 1, 0.5, 'runge-kutta', n=2)
        sol_b_p = self._solve(fb, 2, 3, 1, 0.5, 'trapezoidal')
        self.assertAlmostEqual(sol_b_t, 2.4257869726429155, places=7)
        self.assertAlmostEqual(sol_b_r, 2.481553077697754, places=7)
        self.assertAlmostEqual(sol_b_p, 2.516854718589633, places=7)
        self.assertAlmostEqual(abs(gt_b(3) - sol_b_t), 0.07421302735708446, places=5)
        self.assertAlmostEqual(abs(gt_b(3) - sol_b_r), 0.018446922302246094, places=5)
        self.assertAlmostEqual(abs(gt_b(3) - sol_b_p), 0.016854718589633055, places=5)

        # (c)
        fc = BivariateFunction(lambda x, y: 1 + y / x)
        gt_c = Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2)
        sol_c_t = self._solve(fc, 1, 2, 1, 0.25, 'taylor', n=2)
        sol_c_r = self._solve(fc, 1, 2, 1, 0.25, 'runge-kutta', n=2)
        sol_c_p = self._solve(fc, 1, 2, 1, 0.25, 'trapezoidal')
        self.assertAlmostEqual(sol_c_t, 3.3940488287468473, places=7)
        self.assertAlmostEqual(sol_c_r, 3.372858560090703, places=7)
        self.assertAlmostEqual(sol_c_p, 3.3824397824272534, places=7)
        self.assertAlmostEqual(abs(gt_c(2) - sol_c_t), 1.9922455323730435, places=5)
        self.assertAlmostEqual(abs(gt_c(2) - sol_c_r), 2.0134358010291877, places=5)
        self.assertAlmostEqual(abs(gt_c(2) - sol_c_p), 2.0038545786926374, places=5)


class TestProblemSet8(unittest.TestCase):
    """Tests for Problem Set 8 - Adams-Bashforth, Adams-Moulton, and Predictor-Corrector"""

    def _gt_a(self):
        return (
            Polynomial(0, 1 / 5) * Exponent(Polynomial(0, 3))
            - (1 / 25) * Exponent(Polynomial(0, 3))
            + (1 / 25) * Exponent(Polynomial(0, -2))
        )

    def test_problem1(self):
        """
        Use the two-step Adams-Bashforth explicit method to approximate the
        solution of three IVPs and compute the error at the endpoint.
        (a) y' = t e^{3t} - 2y, y(0) = 0, h = 0.2
        (b) y' = 1 + (t - y)^2, y(2) = 1, h = 0.2
        (c) y' = 1 + y/t, y(1) = 2, h = 0.2
        """
        gt_a = self._gt_a()
        gt_b = Polynomial(0, 1) + 1 / Polynomial(1, -1)
        gt_c = Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2)

        sol_a = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2 * y),
            0, 1, 0
        ).solve(0.2, method='adam-bashforth', step=2, points=[gt_a(0.2)])(1)

        sol_b = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + (t - y) ** 2),
            2, 3, 1
        ).solve(0.2, method='adam-bashforth', step=2, points=[gt_b(2.2)])(3)

        sol_c = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + y / t),
            1, 2, 2
        ).solve(0.2, method='adam-bashforth', step=2, points=[gt_c(1.2)])(2)

        self.assertAlmostEqual(sol_a, 2.8241682513173867, places=7)
        self.assertAlmostEqual(abs(gt_a(1) - sol_a), 0.39493106772210407, places=5)
        self.assertAlmostEqual(sol_b, 2.4884511641911957, places=7)
        self.assertAlmostEqual(abs(gt_b(3) - sol_b), 0.01154883580880428, places=5)
        self.assertAlmostEqual(sol_c, 5.39494164236449, places=7)
        self.assertAlmostEqual(abs(gt_c(2) - sol_c), 0.008647281244599014, places=5)

    def test_problem2(self):
        """
        Redo Problem 1 using the two-step Adams-Moulton implicit method.
        """
        gt_a = self._gt_a()
        gt_b = Polynomial(0, 1) + 1 / Polynomial(1, -1)
        gt_c = Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2)

        sol_a = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2 * y),
            0, 1, 0
        ).solve(0.2, method='adam-moulton', step=2, points=[gt_a(0.2)])(1)

        sol_b = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + (t - y) ** 2),
            2, 3, 1
        ).solve(0.2, method='adam-moulton', step=2, points=[gt_b(2.2)])(3)

        sol_c = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + y / t),
            1, 2, 2
        ).solve(0.2, method='adam-moulton', step=2, points=[gt_c(1.2)])(2)

        self.assertAlmostEqual(sol_a, 3.251286594238682, places=7)
        self.assertAlmostEqual(abs(gt_a(1) - sol_a), 0.0321872751991914, places=5)
        self.assertAlmostEqual(sol_b, 2.499410393690296, places=7)
        self.assertAlmostEqual(abs(gt_b(3) - sol_b), 0.0005896063097039494, places=5)
        self.assertAlmostEqual(sol_c, 5.386529272083939, places=7)
        self.assertAlmostEqual(abs(gt_c(2) - sol_c), 0.00023491096404804068, places=5)

    def test_problem3(self):
        """
        Redo Problem 1 using the three-step Adams-Bashforth and three-step
        Adams-Moulton methods.
        """
        gt_a = self._gt_a()
        gt_b = Polynomial(0, 1) + 1 / Polynomial(1, -1)
        gt_c = Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2)

        # Three-step Adams-Bashforth
        sol_ab_a = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2 * y),
            0, 1, 0
        ).solve(0.2, method='adam-bashforth', step=3, points=[gt_a(0.2), gt_a(0.4)])(1)

        sol_ab_b = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + (t - y) ** 2),
            2, 3, 1
        ).solve(0.2, method='adam-bashforth', step=3, points=[gt_b(2.2), gt_b(2.4)])(3)

        sol_ab_c = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + y / t),
            1, 2, 2
        ).solve(0.2, method='adam-bashforth', step=3, points=[gt_c(1.2), gt_c(1.4)])(2)

        self.assertAlmostEqual(sol_ab_a, 3.0360680209672974, places=7)
        self.assertAlmostEqual(abs(gt_a(1) - sol_ab_a), 0.18303129807219332, places=5)
        self.assertAlmostEqual(sol_ab_b, 2.505134038196207, places=7)
        self.assertAlmostEqual(abs(gt_b(3) - sol_ab_b), 0.005134038196207058, places=5)
        self.assertAlmostEqual(sol_ab_c, 5.384805838441112, places=7)
        self.assertAlmostEqual(abs(gt_c(2) - sol_ab_c), 0.0014885226787786365, places=5)

        # Three-step Adams-Moulton
        sol_am_a = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2 * y),
            0, 1, 0
        ).solve(0.2, method='adam-moulton', step=3, points=[gt_a(0.2), gt_a(0.4)])(1)

        sol_am_b = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + (t - y) ** 2),
            2, 3, 1
        ).solve(0.2, method='adam-moulton', step=3, points=[gt_b(2.2), gt_b(2.4)])(3)

        sol_am_c = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + y / t),
            1, 2, 2
        ).solve(0.2, method='adam-moulton', step=3, points=[gt_c(1.2), gt_c(1.4)])(2)

        self.assertAlmostEqual(sol_am_a, 3.2298092156426033, places=7)
        self.assertAlmostEqual(abs(gt_a(1) - sol_am_a), 0.010709896603112501, places=5)
        self.assertAlmostEqual(sol_am_b, 2.5002011297197355, places=7)
        self.assertAlmostEqual(abs(gt_b(3) - sol_am_b), 0.00020112971973551552, places=5)
        self.assertAlmostEqual(sol_am_c, 5.386255753106247, places=7)
        self.assertAlmostEqual(abs(gt_c(2) - sol_am_c), 3.860801364385935e-05, places=7)

    def test_problem4(self):
        """
        Apply the Adams fourth-order predictor-corrector method with h = 0.2
        (RK4 starting values) to y' = y - t^2 + 1, y(0) = 0.5, 0 <= t <= 2.
        """
        f = BivariateFunction(lambda t, y: y - t ** 2 + 1)
        IVP = FirstOrderLinearODE(f, 0, 2, 0.5)
        sol = IVP.solve(0.2, method='predictor-corrector')
        self.assertAlmostEqual(sol(2), 5.305370671515845, places=7)

    def test_problem5(self):
        """
        Use the Runge-Kutta method of order 2 to approximate the solution of
        the system:
        u1' = 3u1 + 2u2 - (2t^2+1)e^{2t},  u1(0) = 1
        u2' = 4u1 + u2 + (t^2+2t-4)e^{2t}, u2(0) = 0
        on [0, 1] with h = 0.2.
        """
        a, b, h = 0, 1, 0.2
        U0 = Vector(1, 0)
        F = Vector(
            MultiVariableFunction(
                lambda t, u1, u2: 3 * u1 + 2 * u2 - (2 * (t ** 2) + 1) * Exponent(Polynomial(0, 2))(t)
            ),
            MultiVariableFunction(
                lambda t, u1, u2: 4 * u1 + u2 + (t ** 2 + 2 * t - 4) * Exponent(Polynomial(0, 2))(t)
            ),
        )
        IVP = FirstOrderLinearODE(F, a, b, U0)
        sol = IVP.solve(h, method='runge-kutta', n=2)
        self.assertAlmostEqual(sol[-1][0], 3.771001932367687, places=7)
        self.assertAlmostEqual(sol[-1][1], 4.121818316204175, places=7)

    def test_problem6(self):
        """
        Redo Problem 1 using the five-step Adams-Bashforth and five-step
        Adams-Moulton methods (tests the n-step generalization for n > 4).
        (a) y' = t e^{3t} - 2y, y(0) = 0, h = 0.2. Exact: y = (t/5 - 1/25)e^{3t} + (1/25)e^{-2t}.
        (b) y' = 1 + (t - y)^2, y(2) = 1, h = 0.2. Exact: y = t + 1/(1 - t).
        (c) y' = 1 + y/t, y(1) = 2, h = 0.2. Exact: y = t ln(t) + 2t.
        """
        gt_a = self._gt_a()
        gt_b = Polynomial(0, 1) + 1 / Polynomial(1, -1)
        gt_c = Polynomial(0, 1) * Log(Polynomial(0, 1)) + Polynomial(0, 2)

        # Five-step Adams-Bashforth
        sol_ab_a = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2 * y),
            0, 1, 0
        ).solve(0.2, method='adam-bashforth', step=5, points=[gt_a(0.2), gt_a(0.4), gt_a(0.6), gt_a(0.8)])(1)

        sol_ab_b = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + (t - y) ** 2),
            2, 3, 1
        ).solve(0.2, method='adam-bashforth', step=5, points=[gt_b(2.2), gt_b(2.4), gt_b(2.6), gt_b(2.8)])(3)

        sol_ab_c = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + y / t),
            1, 2, 2
        ).solve(0.2, method='adam-bashforth', step=5, points=[gt_c(1.2), gt_c(1.4), gt_c(1.6), gt_c(1.8)])(2)

        self.assertAlmostEqual(sol_ab_a, 3.185400192276253, places=7)
        self.assertAlmostEqual(abs(gt_a(1) - sol_ab_a), 0.03369912676323761, places=5)
        self.assertAlmostEqual(sol_ab_b, 2.501140638997369, places=7)
        self.assertAlmostEqual(abs(gt_b(3) - sol_ab_b), 0.0011406389973691589, places=5)
        self.assertAlmostEqual(sol_ab_c, 5.386217688517386, places=7)
        self.assertAlmostEqual(abs(gt_c(2) - sol_ab_c), 7.667260250521224e-05, places=7)

        # Five-step Adams-Moulton
        sol_am_a = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: t * Exponent(Polynomial(0, 3))(t) - 2 * y),
            0, 1, 0
        ).solve(0.2, method='adam-moulton', step=5, points=[gt_a(0.2), gt_a(0.4), gt_a(0.6), gt_a(0.8)])(1)

        sol_am_b = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + (t - y) ** 2),
            2, 3, 1
        ).solve(0.2, method='adam-moulton', step=5, points=[gt_b(2.2), gt_b(2.4), gt_b(2.6), gt_b(2.8)])(3)

        sol_am_c = FirstOrderLinearODE(
            BivariateFunction(lambda t, y: 1 + y / t),
            1, 2, 2
        ).solve(0.2, method='adam-moulton', step=5, points=[gt_c(1.2), gt_c(1.4), gt_c(1.6), gt_c(1.8)])(2)

        self.assertAlmostEqual(sol_am_a, 3.2202048442288973, places=7)
        self.assertAlmostEqual(abs(gt_a(1) - sol_am_a), 0.001105525189406542, places=5)
        self.assertAlmostEqual(sol_am_b, 2.5000316387876476, places=7)
        self.assertAlmostEqual(abs(gt_b(3) - sol_am_b), 3.1638787647558786e-05, places=7)
        self.assertAlmostEqual(sol_am_c, 5.38629254442049, places=7)
        self.assertAlmostEqual(abs(gt_c(2) - sol_am_c), 1.8166994006918458e-06, places=9)


class TestProblemSet9(unittest.TestCase):
    """Tests for Problem Set 9 - Second-Order BVP: Shooting Methods"""

    def test_problem1(self):
        """
        Apply the Linear Shooting technique with N = 10 to
        y'' = -(2/x)y' + (2/x^2)y + sin(ln x)/x^2, y(1)=1, y(2)=2.
        Exact: y = 1.139x - 0.039/x^2 - 0.3 sin(ln x) - 0.1 cos(ln x).
        """
        a, b, N = 1, 2, 10
        h = (b - a) / (N + 1)
        BVP = SecondOrderLinearODE_BVP(
            Function(lambda x: -(2 / x)),
            Function(lambda x: 2 / x ** 2),
            Function(lambda x: Sin(Log(Polynomial(0, 1)))(x) / (x ** 2)),
            a, b, y0=1, y1=2,
        )
        sol = BVP.solve(h, method='shooting')
        GT = (
            Polynomial(0, 1.139)
            - 0.039 / Polynomial(0, 0, 1)
            - (3 / 10) * Sin(Log(Polynomial(0, 1)))
            - (1 / 10) * Cos(Log(Polynomial(0, 1)))
        )

        expected_sol = [
            1.0, 1.0845389257458282, 1.1702360633395221, 1.2573168688097824,
            1.3458334692760268, 1.4357531316818575, 1.5270029418325044,
            1.6194929848065587, 1.7131285061608457, 1.8078162396379827,
            1.903467576101026,
        ]
        expected_gt = [
            1.0, 1.0840824439316827, 1.169676492566062, 1.2568043588337137,
            1.3454180023626199, 1.4354355137787183, 1.5267607629584172,
            1.6192940433059289, 1.7129377541969864, 1.807599304372986,
            1.9031924463836414,
        ]
        for i in range(N + 1):
            self.assertAlmostEqual(sol(a + i * h), expected_sol[i], places=7)
            self.assertAlmostEqual(GT(a + i * h), expected_gt[i], places=7)

    def test_problem2(self):
        """
        Use the Linear Shooting method with h = 1/4 to solve
        y'' = 4(y - x), y(0) = 0, y(1) = 2.
        Exact: y = e^2 / (e^4 - 1) * (e^{2x} - e^{-2x}) + x.
        """
        a, b, h = 0, 1, 0.25
        BVP = SecondOrderLinearODE_BVP(
            Function(lambda x: 0),
            Function(lambda x: 4),
            Function(lambda x: -4 * x),
            a, b, y0=0, y1=2,
        )
        sol = BVP.solve(h, method='shooting')
        GT = (
            (math.e ** 2 / (math.e ** 4 - 1))
            * (Exponent(Polynomial(0, 2)) - Exponent(Polynomial(0, -2)))
            + Polynomial(0, 1)
        )

        expected_sol = [0.0, 0.45, 0.9, 1.4000000000000001, 2.0]
        expected_gt = [0.0, 0.39367669193066096, 0.8240271368319427, 1.3370861339156979, 2.0]
        n_steps = int((b - a) / h)
        for i in range(n_steps + 1):
            self.assertAlmostEqual(sol(a + i * h), expected_sol[i], places=7)
            self.assertAlmostEqual(GT(a + i * h), expected_gt[i], places=7)

    def test_problem3(self):
        """
        Apply the nonlinear shooting method (with Newton iteration) to
        y'' = (1/8)(32 + 2x^3 - yy'), y(1)=17, y(3)=43/3, N=20.
        Exact: y(x) = x^2 + 16/x.
        """
        a, b, N = 1, 3, 20
        h = (b - a) / (N + 1)
        f = MultiVariableFunction(lambda x, y, z: (1 / 8) * (32 + 2 * x ** 3 - y * z))
        BVP = SecondOrderODE_BVP(f, a, b, y0=17, y1=43 / 3)
        sol = BVP.solve(h, method='shooting_newton', M=100, TOL=1e-5)
        GT = Polynomial(0, 0, 1) + 16 / Polynomial(0, 1)

        expected_sol = [
            17.0, 15.680210825398234, 14.66607059372891, 13.880499376487947,
            13.272192803655221, 12.805505866213215, 12.454809408851439,
            12.201150961314529, 12.030183265415129, 11.930828140110533,
            11.894385968672392, 11.913925363122996, 11.983854601383817,
            12.09961422611917, 12.257452322419587, 12.45425738238914,
            12.687432007771182, 12.954796036696827, 13.254511172977113,
            13.585021531805989, 13.945006108362099,
        ]
        for i in range(N + 1):
            self.assertAlmostEqual(sol(a + i * h), expected_sol[i], places=5)
            self.assertAlmostEqual(GT(a + i * h), expected_sol[i], delta=0.25)

    def test_problem4(self):
        """
        Use the Nonlinear Shooting method with h = 0.5 to approximate
        y'' = -(y')^2 - y + ln x, y(1)=0, y(2)=ln 2.
        Exact: y(x) = ln x.
        """
        a, b, h = 1, 2, 0.5
        f = MultiVariableFunction(lambda x, y, z: -(z ** 2) - y + Log(Polynomial(0, 1))(x))
        BVP = SecondOrderODE_BVP(f, a, b, y0=0, y1=math.log(2))
        sol = BVP.solve(h, method='shooting_newton')
        GT = Log(Polynomial(0, 1))

        expected_sol = [0.0, 0.44606540587021165, 0.6931564654262666]
        expected_gt = [0.0, 0.4054651081081644, 0.6931471805599453]
        n_steps = int((b - a) / h)
        for i in range(n_steps + 1):
            self.assertAlmostEqual(sol(a + i * h), expected_sol[i], places=7)
            self.assertAlmostEqual(GT(a + i * h), expected_gt[i], places=7)


class TestProblemSet10(unittest.TestCase):
    """Tests for Problem Set 10 - Second-Order BVP: Finite Difference Method"""

    def test_problem1(self):
        """
        Use the Linear Finite-Difference Algorithm with N = 9 to approximate
        y'' = -(2/x)y' + (2/x^2)y + sin(ln x)/x^2, y(1)=1, y(2)=2,
        and compare with the linear shooting method.
        """
        a, b, N = 1, 2, 9
        h = (b - a) / (N + 1)
        BVP = SecondOrderLinearODE_BVP(
            Function(lambda x: -2 / x),
            Function(lambda x: 2 / x ** 2),
            Function(lambda x: (Sin(Log(Polynomial(0, 1))) / Polynomial(0, 0, 1))(x)),
            a, b, y0=1, y1=2,
        )
        sol_fd = BVP.solve(h, method='finite_difference')
        sol_sh = BVP.solve(h, method='shooting')

        expected_fd = [
            1.0, 1.092600520720135, 1.1870431287950733, 1.283336870214396,
            1.3814020462286074, 1.4811202621122181, 1.5823598956934564,
            1.6849890183845608, 1.7888817461937458, 1.8939210991571331,
            1.9999999999999993,
        ]
        expected_sh = [
            1.0, 1.0931207035549964, 1.1876172663989932, 1.283787242216369,
            1.3816826919998089, 1.4812447537525195, 1.5823669525635886,
            1.6849259057874906, 1.7887963313670905, 1.8938582193782807,
            2.0,
        ]
        for i in range(N + 2):
            self.assertAlmostEqual(sol_fd(a + i * h), expected_fd[i], places=7)
            self.assertAlmostEqual(sol_sh(a + i * h), expected_sh[i], places=7)

    def test_problem2(self):
        """
        Use the Linear Finite-Difference Algorithm to approximate
        y'' = -(x+1)y' + 2y + (1-x^2)e^{-x}, y(0)=1, y(1)=0, with N=9 and N=19.
        """
        a, b = 0, 1
        N1, N2 = 9, 19
        h1 = (b - a) / (N1 + 1)
        h2 = (b - a) / (N2 + 1)
        BVP = SecondOrderLinearODE_BVP(
            Polynomial(-1, -1),
            Polynomial(2),
            Polynomial(1, 0, -1) * Exponent(Polynomial(0, -1)),
            a, b, y0=1, y1=0,
        )
        sol1 = BVP.solve(h1, method='finite_difference')
        sol2 = BVP.solve(h2, method='finite_difference')

        # Spot-check a few values at the coarser grid spacing h2
        expected_n19 = [
            1.0, 0.8754842088861384, 0.7639151933681421, 0.664215126200481,
            0.5753604389701534, 0.49638248692046416, 0.42636778216583554,
            0.36445782915157593, 0.309848597270125, 0.2617896661758497,
            0.21958307955685577, 0.18258194294380875, 0.1501888005873982,
            0.12185382554752623, 0.09707285594220512, 0.07538530883903964,
            0.056372001575666096, 0.039652908407304364, 0.024884878339682063,
            0.011759337853444506, 1.3924263951508147e-16,
        ]
        for i in range(N2 + 2):
            self.assertAlmostEqual(sol2(a + i * h2), expected_n19[i], places=7)

    def test_problem3(self):
        """
        Use the Linear Finite-Difference Algorithm with N = 4 to approximate
        y'' + 4y = cos x, y(0)=0, y(pi/4)=0.
        Exact: y(x) = -(1/3)cos 2x - (sqrt(2)/6)sin 2x + (1/3)cos x.
        """
        a, b, N = 0, math.pi / 4, 4
        h = (b - a) / (N + 1)
        BVP = SecondOrderLinearODE_BVP(
            Polynomial(0),   # p(x) = 0  (coefficient of y')
            Polynomial(-4),  # q(x) = -4 (coefficient of y, rewriting y''+4y as y''-4y=-cosx)
            Cos(Polynomial(0, 1)),
            a, b, y0=0, y1=0,
        )
        sol = BVP.solve(h, method='finite_difference')
        GT = (
            -(1 / 3) * Cos(Polynomial(0, 2))
            - (math.sqrt(2) / 6) * Sin(Polynomial(0, 2))
            + (1 / 3) * Cos(Polynomial(0, 1))
        )

        expected_approx = [
            0.0, -0.061418449778823124, -0.09240490855208863,
            -0.0908049894567962, -0.05825827234434498, -1.5693826446727677e-17,
        ]
        expected_actual = [
            0.0, -0.06062539597480876, -0.09119580528574917,
            -0.08961337697392985, -0.05749950398798581, -5.551115123125783e-17,
        ]
        for i in range(N + 2):
            self.assertAlmostEqual(sol(a + i * h), expected_approx[i], places=7)
            self.assertAlmostEqual(GT(a + i * h), expected_actual[i], places=7)


class TestProblemSet11(unittest.TestCase):
    """Tests for Problem Set 11 - Linear Systems"""

    def test_problem1(self):
        """
        Use Gaussian elimination with backward substitution to solve
        4x1 - x2 + x3 = 8
        2x1 + 5x2 + 2x3 = 3
        x1 + 2x2 + 4x3 = 11
        Exact solution: x1=1, x2=-1, x3=3.
        """
        A = Matrix(Vector(4, -1, 1), Vector(2, 5, 2), Vector(1, 2, 4))
        b = Vector(8, 3, 11)
        result = LinearSystem(A, b).solve(method='gauss_elimination')
        self.assertAlmostEqual(result[0], 1.0, places=7)
        self.assertAlmostEqual(result[1], -1.0, places=7)
        self.assertAlmostEqual(result[2], 3.0, places=7)

    def test_problem2(self):
        """
        Use Gauss-Jacobi iterative method to solve
        10x1 - x2 + 2x3           =  6
        -x1 + 11x2 - x3 + 3x4    = 25
        2x1 -  x2 + 10x3 - x4    = -11
              3x2 -  x3 + 8x4    = 15
        Unique solution: x = (1, 2, -1, 1).  Stop when ||x^k - x^{k-1}||/||x^k|| < 1e-3.
        """
        A = Matrix(
            Vector(10, -1, 2, 0),
            Vector(-1, 11, -1, 3),
            Vector(2, -1, 10, -1),
            Vector(0, 3, -1, 8),
        )
        b = Vector(6, 25, -11, 15)
        result = LinearSystem(A, b).solve(
            method='gauss_jacobi', TOL=1e-3, initial_approximation=Vector(0, 0, 0, 0)
        )
        self.assertAlmostEqual(result[0], 0.9996741452148707, places=3)
        self.assertAlmostEqual(result[1], 2.0004476715450092, places=3)
        self.assertAlmostEqual(result[2], -1.0003691576845712, places=3)
        self.assertAlmostEqual(result[3], 1.0006191901399695, places=3)

    def test_problem3(self):
        """
        Solve Problem 2 using the Gauss-Seidel iterative method.
        """
        A = Matrix(
            Vector(10, -1, 2, 0),
            Vector(-1, 11, -1, 3),
            Vector(2, -1, 10, -1),
            Vector(0, 3, -1, 8),
        )
        b = Vector(6, 25, -11, 15)
        result = LinearSystem(A, b).solve(
            method='gauss_seidel', TOL=1e-3, initial_approximation=Vector(0, 0, 0, 0)
        )
        self.assertAlmostEqual(result[0], 1.000091280285995, places=3)
        self.assertAlmostEqual(result[1], 2.000021342246459, places=3)
        self.assertAlmostEqual(result[2], -1.0000311471834449, places=3)
        self.assertAlmostEqual(result[3], 0.9999881032596473, places=3)

    def test_problem4(self):
        """
        Use Gauss-Jacobi iterations to attempt solving
        x1 + 2x2 + 3x3 = 5
        2x1 - x2 + 2x3 = 1
        3x1 + x2 - 2x3 = -1
        The method does not converge for this system.
        """
        A = Matrix(Vector(1, 2, 3), Vector(2, -1, 2), Vector(3, 1, -2))
        b = Vector(5, 1, -1)
        result = LinearSystem(A, b).solve(
            method='gauss_jacobi', TOL=1e-3, initial_approximation=Vector(0, 0, 0)
        )
        self.assertIsNone(result)

    def test_problem5(self):
        """
        Use Gauss-Seidel iterations to solve
        7x1 - 2x2 + x3 + 2x4 = 3
        2x1 + 8x2 + 3x3 + x4 = -2
        -x1 + 2x2 + 5x3      = 5
              2x2 - x3 + 4x4 = 4
        """
        A = Matrix(
            Vector(7, -2, 1, 2),
            Vector(2, 8, 3, 1),
            Vector(-1, 2, 5, 0),
            Vector(0, 2, -1, 4),
        )
        b = Vector(3, -2, 5, 4)
        result = LinearSystem(A, b).solve(
            method='gauss_seidel', TOL=1e-3, initial_approximation=Vector(0, 0, 0, 0)
        )
        self.assertAlmostEqual(result[0], -0.46691045573724593, places=3)
        self.assertAlmostEqual(result[1], -0.8083614774349144, places=3)
        self.assertAlmostEqual(result[2], 1.2299624998265166, places=3)
        self.assertAlmostEqual(result[3], 1.7116713636740863, places=3)


class TestProblemSet12(unittest.TestCase):
    """Tests for Problem Set 12 - Nonlinear Finite Difference Method"""

    def test_problem1(self):
        """
        Use the Nonlinear Finite-Difference Algorithm with h = 0.5 to approximate
        y'' = -(y')^2 - y + ln x, y(1) = 0, y(2) = ln 2.
        Exact: y(x) = ln x.
        """
        a, b, h = 1, 2, 0.5
        f = MultiVariableFunction(lambda x, y, z: -(z ** 2) - y + math.log(x))
        BVP = SecondOrderODE_BVP(f, a, b, y0=0, y1=math.log(2))
        sol = BVP.solve(h, method='finite_difference')
        GT = Log(Polynomial(0, 1))

        n_steps = int((b - a) / h)
        for i in range(n_steps + 1):
            self.assertAlmostEqual(sol(a + i * h), GT(a + i * h), delta=1e-2)

    def test_problem2(self):
        """
        Use the Nonlinear Finite-Difference Algorithm with h = 0.25 to approximate
        y'' = -(y')^2 - y + ln x, y(1) = 0, y(2) = ln 2.
        Exact: y(x) = ln x.
        """
        a, b, h = 1, 2, 0.25
        f = MultiVariableFunction(lambda x, y, z: -(z ** 2) - y + math.log(x))
        BVP = SecondOrderODE_BVP(f, a, b, y0=0, y1=math.log(2))
        sol = BVP.solve(h, method='finite_difference')
        GT = Log(Polynomial(0, 1))

        n_steps = int((b - a) / h)
        for i in range(n_steps + 1):
            self.assertAlmostEqual(sol(a + i * h), GT(a + i * h), delta=1e-3)

    def test_problem3(self):
        """
        Apply the Nonlinear Finite-Difference Algorithm to
        y'' = (1/8)(32 + 2x^3 - yy'), y(1) = 17, y(3) = 43/3, N = 20.
        Exact: y(x) = x^2 + 16/x.
        """
        a, b, N = 1, 3, 20
        h = (b - a) / (N + 1)
        f = MultiVariableFunction(lambda x, y, z: (1 / 8) * (32 + 2 * x ** 3 - y * z))
        BVP = SecondOrderODE_BVP(f, a, b, y0=17, y1=43 / 3)
        sol = BVP.solve(h, method='finite_difference', M=100, TOL=1e-5)
        GT = Polynomial(0, 0, 1) + 16 / Polynomial(0, 1)

        for i in range(N + 2):
            self.assertAlmostEqual(sol(a + i * h), GT(a + i * h), delta=1e-2)


if __name__ == '__main__':
    unittest.main()
