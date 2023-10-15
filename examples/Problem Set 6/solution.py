from function import Function, Polynomial, Cos, Sin, Exponent, Log
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
        pass

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
