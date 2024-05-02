from function import Function, Polynomial
import math


class Util:

    @staticmethod
    def choose(s: Function, k: int):
            res = Polynomial(1)
            for i in range(k):
                res *= (-i + s)
            return (1 / math.factorial(k)) * res
    
    @staticmethod
    def delta(k: int, i: int, data: list[tuple[float]]):
        if k == 1:
            return data[i+1][1] - data[i][1]
        return Util.delta(k - 1, i + 1, data) - Util.delta(k - 1, i, data)

    @staticmethod
    def downdelta(k: int, i: int, data: list[tuple[float]]):
        if k == 1:
            return data[i][1] - data[i-1][1]
        return Util.downdelta(k - 1, i, data) - Util.downdelta(k - 1, i - 1, data)
    