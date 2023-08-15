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

    def derivative(self, h=10e-5):
        return Function(lambda x: (self(x + h) - self(x)) / h)
    
    def root(self, *args, **kwargs):
        """
        Wrapper for Solve.root
        """
        from solve import Solve
        return Solve.root(self, *args, **kwargs)
    
    def fixed_point(self, *args, **kwargs):
        """
        Wrapper for Solve.fixed_point
        """
        from solve import Solve
        return Solve.fixed_point(self, *args, **kwargs)

class Polynomial(Function):

    def __init__(self, *coefficients):
        """
            coefficients are in the form a_0, a_1, ... a_n
        """
        self.function = lambda x: sum(a * x ** i for i, a in enumerate(coefficients))
