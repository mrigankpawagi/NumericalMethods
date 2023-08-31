from function import Function, Polynomial
from util import Util

class Solution:

    @staticmethod
    def solve(n):
        return getattr(Solution, f'problem{n}')()

    @staticmethod
    def solve_all(n):
        for i in range(1, n + 1):
            print(f"Problem {i}:")
            print(getattr(Solution, f'problem{i}')())

        # [Solution.solve(i) for i in range(1, n + 1)]

    @staticmethod
    def problem1():
        """Use Newton's forward difference formula to construct interpolating polynomials of degree 1, 2, and 3 for
        the following data"""
        # (a) f(-1/3) if f(-0.75) = -0.07181250, f(-0.5) = -0.02475, f(-0.25) = 0.3349375, f(0) = 1.101000
        data = [(-0.75, -0.07181250), (-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101000)]
        print("Part (a)")
        for n in range(2, 5):
            p1 = Polynomial.interpolate(data[:n], form='forward_diff')
            print(f"n = {n-1}, f(-1/3) = {p1(-1/3)}")

        # (b) f(0.25) if f(0.1) = -0.62049958, f(0.2) = -0.28398668, f(0.3) = 0.00660095, f(0.4) = 0.24842440
        data = [(0.1, -0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)]
        print("\nPart (b)")
        for n in range(2, 5):
            p2 = Polynomial.interpolate(data[:n], form='forward_diff')
            print(f"n = {n-1}, f(0.25) = {p2(0.25)}")
        
        return ""
    
    @staticmethod
    def problem2():
        """Redo problem 1 using Newton's backward difference formula"""
        # (a) f(-1/3) if f(-0.75) = -0.07181250, f(-0.5) = -0.02475, f(-0.25) = 0.3349375, f(0) = 1.101000
        data = [(-0.75, -0.07181250), (-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101000)]
        print("Part (a)")
        for n in range(2, 5):
            p1 = Polynomial.interpolate(data[:n], form='backward_diff')
            print(f"n = {n-1}, f(-1/3) = {p1(-1/3)}")

        # (b) f(0.25) if f(0.1) = -0.62049958, f(0.2) = -0.28398668, f(0.3) = 0.00660095, f(0.4) = 0.24842440
        data = [(0.1, -0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)]
        print("Part (b)")
        for n in range(2, 5):
            p2 = Polynomial.interpolate(data[:n], form='backward_diff')
            print(f"n = {n-1}, f(0.25) = {p2(0.25)}")
        
        return ""

    @staticmethod
    def problem3():
        """Find the degree of the polynomial which interpolates the following data.
        f(-2) = 1, f(-1) = 4, f(0) = 11, f(1) = 16, f(2) = 13, f(3) = -4"""
        data = [(-2, 1), (-1, 4), (0, 11), (1, 16), (2, 13), (3, -4)]
        p = Polynomial.interpolate(data)
        print(f"Forward differences: {Util.delta(1, 0, data)}, {Util.delta(2, 0, data)}, {Util.delta(3, 0, data)}, {Util.delta(4, 0, data)}, {Util.delta(5, 0, data)}")
        return "The degree of the polynomial is 3\n"

    @staticmethod
    def problem4():
        """
        Use appropriate Lagrange interpolating polynomials of degree 1, 2, and 3 to
        approximate the data given in Problem 1."""
        # (a) f(-1/3) if f(-0.75) = -0.07181250, f(-0.5) = -0.02475, f(-0.25) = 0.3349375, f(0) = 1.101000
        data = [(-0.75, -0.07181250), (-0.5, -0.02475), (-0.25, 0.3349375), (0, 1.101000)]
        print("Part (a)")
        for n in range(2, 5):
            p1 = Polynomial.interpolate(data[:n], method='lagrange')
            print(f"n = {n-1}, f(-1/3) = {p1(-1/3)}")

        # (b) f(0.25) if f(0.1) = -0.62049958, f(0.2) = -0.28398668, f(0.3) = 0.00660095, f(0.4) = 0.24842440
        data = [(0.1, -0.62049958), (0.2, -0.28398668), (0.3, 0.00660095), (0.4, 0.24842440)]
        print("\nPart (b)")
        for n in range(2, 5):
            p2 = Polynomial.interpolate(data[:n], method='lagrange')
            print(f"n = {n-1}, f(0.25) = {p2(0.25)}")
        
        return ""
