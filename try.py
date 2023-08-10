from function import Function
from solve import Solve

f = Function.Polynomial(-0.3, 0.4, 1)

print(Solve.root(f, method="bisection", a = 0.2, b = 0.6))
print(Solve.root(f, method="fixed_point", p0 = 0.4))
