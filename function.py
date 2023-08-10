class Function:

    @staticmethod
    def Polynomial(*coefficients):
        """
            coefficients are in the form a_0, a_1, ... a_n
        """
        return lambda x: sum(a * x ** i for i, a in enumerate(coefficients))

    @staticmethod
    def FixedPoint(f):
        return lambda x: x - f(x)
