# Numerical Methods

[![Tests](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml/badge.svg)](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)

A Python library implementing a broad set of classical numerical methods — from root finding and interpolation to ODE solvers and linear systems — built around a clean, composable function abstraction.

Developed as part of UMC 202 at the Indian Institute of Science, Bengaluru.

## Features

### Functions and Calculus
- **Function arithmetic** — add, subtract, multiply, divide, compose, and exponentiate functions
- **Built-in types** — `Polynomial`, `Exponent`, `Sin`, `Cos`, `Tan`, `Log`
- **Differentiation** — forward, backward, and central differences; *n*-th order derivatives
- **Integration** — rectangular, midpoint, trapezoidal, Simpson's, and Gaussian quadrature (*n* = 1, 2)
- **Multivariable functions** — `MultiVariableFunction`, `BivariateFunction`
- **Linear algebra** — `Vector` and `Matrix` types with standard arithmetic

### Root Finding
| Method | Class/Call |
|---|---|
| Bisection | `f.root('bisection', a=..., b=...)` |
| Newton-Raphson | `f.root('newton', p0=...)` |
| Modified Newton | `f.root('modified_newton', p0=...)` |
| Secant | `f.root('secant', p0=..., p1=...)` |
| Regula Falsi | `f.root('regula_falsi', p0=..., p1=...)` |
| Fixed-Point Iteration | `f.fixed_point(p0=...)` |

### Interpolation
`Polynomial.interpolate(data, method=...)` supports:
- Lagrange interpolation
- Newton's divided difference
- Newton forward / backward difference

### ODE Solvers
| Problem Type | Class | Methods |
|---|---|---|
| 1st-order IVP: *y′ = f(x, y)* | `FirstOrderLinearODE` | Euler, Taylor (*n*=1,2), Runge-Kutta (*n*=1–4), Trapezoidal, Adams-Bashforth (*n*=2–4), Adams-Moulton (*n*=2–4), Predictor-Corrector |
| 2nd-order linear BVP | `SecondOrderLinearODE_BVP` | Shooting, Finite Difference |
| 2nd-order nonlinear IVP | `SecondOrderODE_IVP` | All 1st-order methods |
| 2nd-order nonlinear BVP | `SecondOrderODE_BVP` | Shooting (Newton iteration) |

### Linear Systems
`LinearSystem(A, b).solve(method=...)` supports:
- Gaussian elimination with backward substitution
- Gauss-Jacobi iteration
- Gauss-Seidel iteration

## Usage

```python
from numericalmethods import Polynomial, Sin, FirstOrderLinearODE, BivariateFunction, LinearSystem, Vector, Matrix

# Root finding
f = Polynomial(-6, 14, -7, 1)          # x^3 - 7x^2 + 14x - 6
root = f.root('bisection', a=0, b=1, TOLERANCE=1e-6)

# Integration
g = Sin(Polynomial(0, 1))              # sin(x)
area = g.integrate(0, 3.14159, method='simpson', n=100)

# Polynomial interpolation
data = [(0, 1), (1, 2.7183), (2, 7.389), (3, 20.086)]
p = Polynomial.interpolate(data, method='newton')

# Solve y' = y, y(0) = 1 on [0, 1]
f = BivariateFunction(lambda x, y: y)
ode = FirstOrderLinearODE(f, a=0, b=1, y0=1)
approx = ode.solve(h=0.1, method='runge-kutta', n=4)

# Linear system Ax = b
A = Matrix(Vector(2, 1), Vector(5, 7))
b = Vector(11, 13)
x = LinearSystem(A, b).solve()
```

Plotting requires `matplotlib` and is available on any `Function` via `f.plot(min, max)`.

## Running Tests

```bash
python -m unittest discover -s tests -v
```
