from function import Matrix, Vector, LinearSystem
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
        Use Gaussian elimination with backward substitution with tolerance 10^-2 to
        solve the following linear system
        4x1 - x2 + x3 = 8
        2x1 + 5x2 + 2x3 = 3
        x1 + 2x2 + 4x3 = 11
        The exact solution of the system is x1 = 1, x2 = -1, x3 = 3.
        """
        A = Matrix(
            Vector(4, -1, 1),
            Vector(2, 5, 2),
            Vector(1, 2, 4)
        )
        b = Vector(8, 3, 11)
        system = LinearSystem(A, b)
        return system.solve(method='gauss_elimination')
    
    @staticmethod
    def problem2():
        """
        The following linear system
        10x1 - x2 + 2x3 = 6
        -x1 + 11x2 - x3 + 3x4 = 25
        2x1 - x2 + 10x3 - x4 = -11
        3x2 - x3 + 8x4 = 15
        has the unique solution x = (1, 2, -1, 1). Use the Gauss Jacobi's iterative technique
        to find the approximations x^(k) to x with x^(0) = (0, 0, 0, 0) until
        ||x^(k) - x^(k-1)||/||x^(k)|| < 10^-3
        where ||x|| = max{|x1|, |x2|, ..., |x4|}.
        """
        A = Matrix(
            Vector(10, -1, 2, 0),
            Vector(-1, 11, -1, 3),
            Vector(2, -1, 10, -1),
            Vector(0, 3, -1, 8)
        )
        b = Vector(6, 25, -11, 15)
        system = LinearSystem(A, b)
        return system.solve(method='gauss_jacobi', TOL=1e-3, initial_approximation=Vector(0, 0, 0, 0))
    
    @staticmethod
    def problem3():
        """
        Solve problem 2 by Gauss Seidel iterative technique. 
        """
        A = Matrix(
            Vector(10, -1, 2, 0),
            Vector(-1, 11, -1, 3),
            Vector(2, -1, 10, -1),
            Vector(0, 3, -1, 8)
        )
        b = Vector(6, 25, -11, 15)
        system = LinearSystem(A, b)
        return system.solve(method='gauss_seidel', TOL=1e-3, initial_approximation=Vector(0, 0, 0, 0))

    @staticmethod
    def problem4():
        """
        Use Gauss-Jacobi Iterations to attempt solving the linear system
        x1 + 2x2 + 3x3 = 5
        2x1 - x2 + 2x3 = 1
        3x1 + x2 - 2x3 = -1
        """
        A = Matrix(
            Vector(1, 2, 3),
            Vector(2, -1, 2),
            Vector(3, 1, -2)
        )
        b = Vector(5, 1, -1)
        system = LinearSystem(A, b)
        sol = system.solve(method='gauss_jacobi', TOL=1e-3, initial_approximation=Vector(0, 0, 0))
        if sol is None:
            return "Could not converge."
        return sol
    
    @staticmethod
    def problem5():
        """
        Use Gauss-Seidel Iterations to attempt solving the linear system
        2x1 + 8x2 + 3x3 + x4 = -2
        2x2 - x3 + 4x4 = 4
        7x1 - 2x2 + x3 + 2x4 = 3
        -x1 + 2x2 + 5x3 = 5
        """
        A = Matrix(
            Vector(7, -2, 1, 2),
            Vector(2, 8, 3, 1),
            Vector(-1, 2, 5, 0),
            Vector(0, 2, -1, 4)
        )
        b = Vector(3, -2, 5, 4)
        system = LinearSystem(A, b)
        sol = system.solve(method='gauss_seidel', TOL=1e-3, initial_approximation=Vector(0, 0, 0, 0))
        if sol is None:
            return "Could not converge."
        return sol
