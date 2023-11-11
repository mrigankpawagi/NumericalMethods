from function import Function, Polynomial, Cos, Sin, Exponent, Log, SecondOrderLinearODE_BVP, MultiVariableFunction, SecondOrderODE_BVP
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
        Apply the Linear Shooting technique with N = 10 to the boundary value problem
        y'' = -(2/x)y' + (2/x^2)y + sin(ln(x))/x^2, 1 <= x <= 2, y(1) = 1, y(2) = 2 
        and compare the result to those of the exact solution
        y = c1x + c2/x^2 -(3/10)sin(ln(x)) - (1/10)cos(ln(x)) where c1 = 1.139 and c2 = -0.039
        """
        a = 1
        b = 2
        N = 10
        h = (b - a) / (N+1)
        
        BVP = SecondOrderLinearODE_BVP(
            Function(lambda x: -(2/x)),
            Function(lambda x: (2/x**2)),
            Function(lambda x: Sin(Log(Polynomial(0, 1)))(x) / (x**2)),
            a, b, y0=1, y1=2
        )
        
        sol = BVP.solve(h, method='shooting')
        
        GT = Polynomial(0, 1.139) - 0.039 / Polynomial(0, 0, 1) - (3/10) * Sin(Log(Polynomial(0, 1))) - (1/10) * Cos(Log(Polynomial(0, 1)))

        return {
            'solution': [sol(a + i * h) for i in range(N + 1)],
            'ground_truth': [GT(a + i * h) for i in range(N + 1)],
            'errors': [abs(sol(a + i * h) - GT(a + i * h)) for i in range(N + 1)]
        }
        
        
    
    @staticmethod
    def problem2():
        """
        The boundary value problem 
        y'' = 4(y-x), 0 <= x <= 1, y(0) = 0, y(1) = 2
        has the solution y(x) = (e^2)(e^4 - 1)^(-1)(e^(2x) - e^(-2x)) + x. Use the Linear
        Shooting method to approximate the solution and compare the result to the actual
        solution for h = 1/4.
        """
        a = 0
        b = 1
        h = 1/4
        
        BVP = SecondOrderLinearODE_BVP(
            Function(lambda x: 0),
            Function(lambda x: 4),
            Function(lambda x: -4 * x),
            a, b, y0=0, y1=2
        )
        
        sol = BVP.solve(h, method='shooting')
        
        GT = ((math.e ** 2) / (math.e ** 4 - 1)) * (Exponent(Polynomial(0, 2)) - Exponent(Polynomial(0, -2))) + Polynomial(0, 1)
        
        return {
            'solution': [sol(a + i * h) for i in range(int((b - a) / h) + 1)],
            'ground_truth': [GT(a + i * h) for i in range(int((b - a) / h) + 1)],
            'errors': [abs(sol(a + i * h) - GT(a + i * h)) for i in range(int((b - a) / h) + 1)]
        }
    
    @staticmethod
    def problem3():
        """
        Apply the shooting method with Newton's method to the boundary value problem
        y'' = (1/8)(32 + 2x^3 - yy'), 1 <= x <= 3, y(1) = 17, y(3) = 43/3.
        Use N = 20, M = 10 and TOL = 10^(-5), and compare the results with the exact solution
        y(x) = x^2 + 16/x.
        """
        a = 1
        b = 3
        N = 20
        h = (b - a) / (N + 1)
             
        f = MultiVariableFunction(lambda x, y, z: (1/8) * (32 + 2 * x**3 - y * z))
        BVP = SecondOrderODE_BVP(f, a, b, y0=17, y1=43/3)
        sol = BVP.solve(h, method='shooting_newton', M=100, TOL=1e-5)
        
        GT = Polynomial(0, 0, 1) + 16 / Polynomial(0, 1)

        return {
            'solution': [sol(a + i * h) for i in range(N + 1)],
            'ground_truth': [GT(a + i * h) for i in range(N + 1)],
            'errors': [abs(sol(a + i * h) - GT(a + i * h)) for i in range(N + 1)]
        }
    
    @staticmethod
    def problem4():
        """
        Use the Nonlinear Shooting method with h = 0.5 to approximate the solution to 
        the boundary value problem
        y'' = -(y')**2 - y + lnx, 1 <= x <= 2, y(1) = 0, y(2) = ln2
        
        Compare the results to the actual solution y(x) = ln(x).
        """
        a = 1
        b = 2
        h = 0.5
        
        f = MultiVariableFunction(lambda x, y, z: -(z**2) - y + Log(Polynomial(0, 1))(x))
        BVP = SecondOrderODE_BVP(f, a, b, y0=0, y1=math.log(2))
        sol = BVP.solve(h, method='shooting_newton')
        
        GT = Log(Polynomial(0, 1))
        
        return {
            'solution': [sol(a + i * h) for i in range(int((b - a) / h) + 1)],
            'ground_truth': [GT(a + i * h) for i in range(int((b - a) / h) + 1)],
            'errors': [abs(sol(a + i * h) - GT(a + i * h)) for i in range(int((b - a) / h) + 1)]
        }
