from function import Function, Sin, Cos, Log, Polynomial, Exponent, SecondOrderLinearODE_BVP
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
        Use the Linear Finite-Difference Algorithm with N = 9 to approximate
        the solution to the boundary value problem
        y'' = -(2/x)y' + (2/x^2)y + sin(lnx)/x^2, 1 <= x <= 2, y(1) = 1, y(2) = 2
        and compare the results to those obtained using the linear shooting method for 
        the same problem. 
        """
        a = 1
        b = 2
        N = 9
        h = (b - a)/(N + 1)
        
        BVP = SecondOrderLinearODE_BVP(
            Function(lambda x: -2/x),
            Function(lambda x: 2/x**2),
            Function(lambda x: (Sin(Log(Polynomial(0, 1)))/Polynomial(0, 0, 1))(x)),
            a, b, y0=1, y1=2
        )
        sol1 = BVP.solve(h, method='shooting')
        sol2 = BVP.solve(h, method='finite_difference')
        
        return {
            "finite_difference": [sol2(a + i * h) for i in range(N + 2)],
            "shooting": [sol1(a + i * h) for i in range(N + 2)],
            "difference": [abs(sol2(a + i * h) - sol1(a + i * h)) for i in range(N + 2)]
        }
        
    
    @staticmethod
    def problem2():
        """
        Consider the boundary value problem
        y'' = -(x + 1)y' + 2y + (1 - x^2)e^(-x), 0 <= x <= 1, y(0) = 1, y(1) = 0.
        Use N = 9 and N = 19 respectively and apply the Linear Finite-Difference Algorithm
        to approximate the solution to the above boundary value problem.
        """
        a = 0
        b = 1
        N1 = 9
        N2 = 19
        h1 = (b - a)/(N1 + 1)
        h2 = (b - a)/(N2 + 1)
        
        BVP = SecondOrderLinearODE_BVP(
            Polynomial(-1, -1),
            Polynomial(2),
            Polynomial(1, 0, -1) * Exponent(Polynomial(0, -1)),
            a, b, y0=1, y1=0 
        )
        sol1 = BVP.solve(h1, method='finite_difference')
        sol2 = BVP.solve(h2, method='finite_difference')
        
        return {
            "N=9": [sol1(a + i * h2) for i in range(N2 + 2)],
            "N=19": [sol2(a + i * h2) for i in range(N2 + 2)],
            "difference": [abs(sol1(a + i * h2) - sol2(a + i * h2)) for i in range(N2 + 2)]
        }
    
    @staticmethod
    def problem3():
        """
        Use the Linear Finite-Difference Algorithm with N = 4 to approximate the solution of the
        boundary value problem,
        y'' + 4y = cos x, 0 <= x <= pi/4, y(0) = 0, y(pi/4) = 0,
        and compare the results to the actual solution where the actual solution is given by
        y(x) = -1/3 cos 2x - sqrt(2)/6 sin 2x + 1/3 cos x.
        """
        a = 0
        b = math.pi/4
        N = 4
        h = (b - a)/(N + 1)
        
        BVP = SecondOrderLinearODE_BVP(
            Polynomial(0),
            Polynomial(-4),
            Cos(Polynomial(0, 1)),
            a, b, y0=0, y1=0
        )
        
        sol = BVP.solve(h, method='finite_difference')
        
        GT = - (1/3) * Cos(Polynomial(0, 2)) - (math.sqrt(2)/6) * Sin(Polynomial(0, 2)) + (1/3) * Cos(Polynomial(0, 1))
        
        return {
            "approximation": [sol(a + i * h) for i in range(N + 2)],
            "actual": [GT(a + i * h) for i in range(N + 2)],
            "difference": [abs(sol(a + i * h) - GT(a + i * h)) for i in range(N + 2)]
        }
