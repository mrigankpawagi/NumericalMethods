# Numerical Methods

[![Tests](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml/badge.svg)](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)

A Python library implementing a broad set of classical numerical methods — from root finding and interpolation to ODE solvers and linear systems — built around a clean, composable function abstraction.

Developed as part of UMC 202 at the Indian Institute of Science, Bengaluru.

## Package Layout

```
numericalmethods/
├── enums.py              # Typed enums for all method selectors
├── functions/
│   ├── base.py           # Function — core callable abstraction
│   ├── elementary.py     # Polynomial, Exponent, Sin, Cos, Tan, Log
│   └── multivariate.py   # MultiVariableFunction, BivariateFunction
├── linalg/
│   ├── vector.py         # Vector
│   ├── matrix.py         # Matrix
│   └── linear_system.py  # LinearSystem
├── ode/
│   ├── base.py           # OrdinaryDifferentialEquation, LinearODE
│   ├── first_order.py    # FirstOrderLinearODE
│   └── second_order.py   # SecondOrderLinearODE_BVP, SecondOrderODE_IVP/BVP
└── util.py               # Util (forward/backward differences, binomial helper)
```

Everything is re-exported from the top-level `numericalmethods` package, so a single import line is all you need.

## Features

### Functions and Calculus
- **Function arithmetic** — add, subtract, multiply, divide, compose, and exponentiate functions
- **Built-in types** — `Polynomial`, `Exponent`, `Sin`, `Cos`, `Tan`, `Log`
- **Differentiation** — forward, backward, and central differences; *n*-th order derivatives
- **Integration** — rectangular, midpoint, trapezoidal, Simpson's, and Gaussian quadrature (*n* = 1, 2)
- **Multivariable functions** — `MultiVariableFunction`, `BivariateFunction`
- **Linear algebra** — `Vector` and `Matrix` types with standard arithmetic

### Enums

All method selectors are typed enums that also compare equal to their string equivalents, preserving backward compatibility.

| Enum | Values |
|---|---|
| `DifferentiationMethod` | `FORWARD`, `BACKWARD`, `CENTRAL` |
| `IntegrationMethod` | `RECTANGULAR`, `MIDPOINT`, `TRAPEZOIDAL`, `SIMPSON`, `GAUSS` |
| `RootFindingMethod` | `BISECTION`, `NEWTON`, `MODIFIED_NEWTON`, `SECANT`, `REGULA_FALSI` |
| `InterpolationMethod` | `LAGRANGE`, `NEWTON` |
| `InterpolationForm` | `STANDARD`, `FORWARD_DIFF`, `BACKWARD_DIFF` |
| `ODEMethod` | `EULER`, `RUNGE_KUTTA`, `TAYLOR`, `TRAPEZOIDAL`, `ADAMS_BASHFORTH`, `ADAMS_MOULTON`, `PREDICTOR_CORRECTOR` |
| `BVPMethod` | `SHOOTING`, `FINITE_DIFFERENCE` |
| `NonlinearBVPMethod` | `SHOOTING_NEWTON`, `FINITE_DIFFERENCE` |
| `LinearSolverMethod` | `GAUSS_ELIMINATION`, `GAUSS_JACOBI`, `GAUSS_SEIDEL` |

### Root Finding
| Method | Enum | String alias |
|---|---|---|
| Bisection | `RootFindingMethod.BISECTION` | `"bisection"` |
| Newton-Raphson | `RootFindingMethod.NEWTON` | `"newton"` |
| Modified Newton | `RootFindingMethod.MODIFIED_NEWTON` | `"modified_newton"` |
| Secant | `RootFindingMethod.SECANT` | `"secant"` |
| Regula Falsi | `RootFindingMethod.REGULA_FALSI` | `"regula_falsi"` |
| Fixed-Point Iteration | — | `f.fixed_point(p0=...)` |

### Interpolation
`Polynomial.interpolate(data, method=..., form=...)` supports:
- Lagrange interpolation (`InterpolationMethod.LAGRANGE`)
- Newton's divided difference (`InterpolationMethod.NEWTON`, `InterpolationForm.STANDARD`)
- Newton forward difference (`InterpolationMethod.NEWTON`, `InterpolationForm.FORWARD_DIFF`)
- Newton backward difference (`InterpolationMethod.NEWTON`, `InterpolationForm.BACKWARD_DIFF`)

### ODE Solvers
| Problem Type | Class | Methods |
|---|---|---|
| 1st-order IVP: *y′ = f(x, y)* | `FirstOrderLinearODE` | `EULER`, `TAYLOR` (*n*=1,2), `RUNGE_KUTTA` (*n*=1–4), `TRAPEZOIDAL`, `ADAMS_BASHFORTH`, `ADAMS_MOULTON`, `PREDICTOR_CORRECTOR` |
| 2nd-order linear BVP | `SecondOrderLinearODE_BVP` | `SHOOTING`, `FINITE_DIFFERENCE` |
| 2nd-order nonlinear IVP | `SecondOrderODE_IVP` | All 1st-order methods |
| 2nd-order nonlinear BVP | `SecondOrderODE_BVP` | `SHOOTING_NEWTON`, `FINITE_DIFFERENCE` |

### Linear Systems
`LinearSystem(A, b).solve(method=...)` supports:
- `LinearSolverMethod.GAUSS_ELIMINATION` — Gaussian elimination with backward substitution
- `LinearSolverMethod.GAUSS_JACOBI` — Gauss-Jacobi iteration
- `LinearSolverMethod.GAUSS_SEIDEL` — Gauss-Seidel iteration

## Usage

```python
from numericalmethods import (
    Polynomial, Sin, BivariateFunction,
    FirstOrderLinearODE, LinearSystem, Vector, Matrix,
    RootFindingMethod, IntegrationMethod, ODEMethod, LinearSolverMethod,
)

# Root finding — enum or string both work
f = Polynomial(-6, 14, -7, 1)          # x³ - 7x² + 14x - 6
root = f.root(RootFindingMethod.BISECTION, a=0, b=1, TOLERANCE=1e-6)

# Integration
g = Sin(Polynomial(0, 1))              # sin(x)
area = g.integrate(0, 3.14159, method=IntegrationMethod.SIMPSON, n=100)

# Polynomial interpolation
data = [(0, 1), (1, 2.7183), (2, 7.389), (3, 20.086)]
p = Polynomial.interpolate(data)       # Newton divided-difference by default

# Solve y' = y, y(0) = 1 on [0, 1]
f = BivariateFunction(lambda x, y: y)
ode = FirstOrderLinearODE(f, a=0, b=1, y0=1)
approx = ode.solve(h=0.1, method=ODEMethod.RUNGE_KUTTA, n=4)

# Linear system Ax = b
A = Matrix(Vector(2, 1), Vector(5, 7))
b = Vector(11, 13)
x = LinearSystem(A, b).solve(LinearSolverMethod.GAUSS_ELIMINATION)
```

Plotting requires `matplotlib` and is available on any `Function` via `f.plot(min, max)`.

## Running Tests

```bash
python -m unittest discover -s tests -v
```

