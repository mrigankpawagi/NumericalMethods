from function import Exponent, Polynomial

# Find global minimum of e^x - x
f = Exponent(Polynomial(0, 1)) - Polynomial(0, 1)
g = -1 + Exponent(Polynomial(0, 1)) # derivative of e^x - x
print(g.root(method='newton', p0=0.5, return_iterations=True))

g.differentiate(Exponent(Polynomial(0, 1)))
print(g.root(method='newton', p0=0.5, return_iterations=True))
# print(f(g.root(method='newton', p0=0.5)))

# Find global minimum of e^x - x^3
f = Exponent(Polynomial(0, 1)) - Polynomial(0, 0, 0, 1)
g = Exponent(Polynomial(0, 1)) - Polynomial(0, 0, 3) # derivative of e^x - x^3
print(g.root(method='newton', p0=4, return_iterations=True))
print(f(g.root(method='newton', p0=4)))

g.differentiate(Exponent(Polynomial(0, 1)) - Polynomial(0, 6))
# print(g.root(method='newton', p0=4, return_iterations=True))

# Find P_n to Interpolate f(x) = 1/(1+x^2) with x in [-1, 1] at equally spaced points for different values of n
f = 1 / Polynomial(1, 0, 1)
min = -10
max = 10

for n in range(2, 15):
    # get n equally spaced points in [min, max] including the end points
    x = [min + i * (max - min) / (n - 1) for i in range(n)]
    P = Polynomial.interpolate(x, f=f)
    f.plot(min=-10, max=10, clear=True)
    P.plot(min=-10, max=10, file=f"fig{n}.png")
