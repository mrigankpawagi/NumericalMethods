class Function:

    def __init__(self, function):
        self.function = function

    def __call__(self, x):
        return self.function(x)

    def __truediv__(f, g):
        return Function(lambda x: f(x) / g(x))
    
    def __mul__(f, g):
        return Function(lambda x: f(x) * g(x))
    
    def __add__(f, g):
        return Function(lambda x: f(x) + g(x))
    
    def __sub__(f, g):
        return Function(lambda x: f(x) - g(x))

class Polynomial(Function):

        def __init__(self, *coefficients):
            """
                coefficients are in the form a_0, a_1, ... a_n
            """
            self.function = lambda x: sum(a * x ** i for i, a in enumerate(coefficients))
