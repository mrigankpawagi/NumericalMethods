# Documentation

## Table of Contents

- [Installation](#installation)
- [Package layout](#package-layout)
- [Functions](#functions)
  - [Arithmetic composition](#arithmetic-composition)
  - [Differentiation](#differentiation)
  - [Integration](#integration)
  - [Root finding](#root-finding)
  - [Fixed-point iteration](#fixed-point-iteration)
  - [Interpolation](#interpolation)
- [Linear algebra](#linear-algebra)
  - [Vector and Matrix](#vector-and-matrix)
  - [Solving linear systems](#solving-linear-systems)
- [Ordinary differential equations](#ordinary-differential-equations)
  - [First-order IVP](#first-order-ivp)
  - [Second-order BVP (linear)](#second-order-bvp-linear)
  - [Second-order BVP (nonlinear)](#second-order-bvp-nonlinear)
  - [Second-order IVP](#second-order-ivp)
- [Utilities](#utilities)
- [Method reference](#method-reference)

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/mrigankpawagi/PySolverKit.git
cd PySolverKit
pip install -e .
```

---

## Package layout

```
pysolverkit/
├── enums.py              # Method selectors
├── functions/
│   ├── base.py           # Function abstraction
│   ├── elementary.py     # Polynomial, Exponent, Sin, Cos, Tan, Log
│   └── multivariate.py   # MultiVariableFunction, BivariateFunction
├── linalg/
│   ├── vector.py         # Vector
│   ├── matrix.py         # Matrix
│   └── linear_system.py  # LinearSystem
├── ode/
│   ├── base.py           # ODE base classes
│   ├── first_order.py    # FirstOrderLinearODE
│   └── second_order.py   # Second-order IVP/BVP solvers
└── util.py               # Utility helpers
```

Everything is re-exported from the top-level `pysolverkit` package.

---

## Functions

All function types inherit from `Function` and support arithmetic composition, differentiation, integration, and root finding.

### Built-in function types

| Class | Description | Example |
|---|---|---|
| `Function(f)` | Wraps any callable `f(x) -> float` | `Function(math.sin)` |
| `Polynomial(a0, a1, ..., an)` | Polynomial with coefficients in ascending degree | `Polynomial(-6, 14, -7, 1)` → x³ − 7x² + 14x − 6 |
| `Exponent(f)` | e^f(x) | `Exponent(Polynomial(0, -1))` → e^(−x) |
| `Sin(f)` | sin(f(x)) | `Sin(Polynomial(0, 1))` → sin(x) |
| `Cos(f)` | cos(f(x)) | `Cos(Polynomial(0, 2))` → cos(2x) |
| `Tan(f)` | tan(f(x)) | `Tan(Polynomial(0, 1))` |
| `Log(f)` | ln(f(x)) | `Log(Polynomial(0, 1))` → ln(x) |

### Arithmetic composition

Functions compose naturally with operators:

```python
from pysolverkit import Polynomial, Exponent

f = Polynomial(0, 1)         # x
g = Exponent(Polynomial(0, -1))  # e^(-x)

h = f - g                    # x - e^(-x)
k = f * g                    # x * e^(-x)
```

### Differentiation

```python
from pysolverkit import Polynomial, DifferentiationMethod

f = Polynomial(0, 0, 1)  # x^2
df = f.differentiate()   # Forward difference (default), step h=1e-5

# Choose scheme
df_central = f.differentiate(method=DifferentiationMethod.CENTRAL)

# Set exact derivative
f.differentiate(Polynomial(0, 2))  # stores 2x as the exact derivative
```

Available methods: `DifferentiationMethod.FORWARD`, `BACKWARD`, `CENTRAL`.

### Integration

```python
from pysolverkit import Sin, Polynomial, IntegrationMethod
import math

f = Sin(Polynomial(0, 1))  # sin(x)

# Numerical quadrature
area = f.integrate(0, math.pi, method=IntegrationMethod.SIMPSON, n=100)

# Set and use exact antiderivative
f.integrate(lambda x: -math.cos(x))          # stores exact antiderivative
exact_area = f.integrate(0, math.pi)          # uses exact antiderivative
```

Available methods: `IntegrationMethod.RECTANGULAR`, `MIDPOINT`, `TRAPEZOIDAL`, `SIMPSON`, `GAUSS`.

### Root finding

```python
from pysolverkit import Polynomial, Exponent, RootFindingMethod

f = Polynomial(-6, 14, -7, 1)  # x^3 - 7x^2 + 14x - 6

# Bisection — bracket [a, b]
root = f.root(RootFindingMethod.BISECTION, a=0, b=1, TOLERANCE=1e-6)

# Newton's method — initial guess p0 (requires derivative set first)
f.differentiate(Polynomial(14, -14, 3))
root = f.root(RootFindingMethod.NEWTON, p0=0.5, TOLERANCE=1e-6)

# Secant method — two initial guesses
root = f.root(RootFindingMethod.SECANT, p0=0.0, p1=1.0, TOLERANCE=1e-6)

# Regula Falsi — bracket [p0, p1]
root = f.root(RootFindingMethod.REGULA_FALSI, p0=0.0, p1=1.0, TOLERANCE=1e-6)

# Modified Newton — repeated root
root = f.root(RootFindingMethod.MODIFIED_NEWTON, p0=0.5, TOLERANCE=1e-6)

# Return iteration count alongside result
root, iters = f.root(RootFindingMethod.BISECTION, a=0, b=1, TOLERANCE=1e-4, return_iterations=True)
```

Available methods: `RootFindingMethod.BISECTION`, `NEWTON`, `MODIFIED_NEWTON`, `SECANT`, `REGULA_FALSI`.

### Fixed-point iteration

```python
from pysolverkit import Function

g = Function(lambda x: (x + 2/x) / 2)   # Babylonian sqrt(2) iteration
result = g.fixed_point(p0=1.0, TOLERANCE=1e-8)
```

### Interpolation

```python
from pysolverkit import Polynomial, InterpolationMethod, InterpolationForm

# Lagrange interpolation from data points
xs = [0.0, 0.5, 1.0]
ys = [1.0, 1.6487, 2.7183]
p = Polynomial.interpolate(xs, ys, method=InterpolationMethod.LAGRANGE)

# Newton's divided-difference form
p = Polynomial.interpolate(xs, ys, method=InterpolationMethod.NEWTON)

# Newton's forward/backward difference form (equal spacing required)
p = Polynomial.interpolate(xs, ys, method=InterpolationMethod.NEWTON,
                           form=InterpolationForm.FORWARD_DIFF)
```

Available methods: `InterpolationMethod.LAGRANGE`, `NEWTON`.  
Available forms: `InterpolationForm.STANDARD`, `FORWARD_DIFF`, `BACKWARD_DIFF`.

---

## Linear algebra

### Vector and Matrix

```python
from pysolverkit import Vector, Matrix

v = Vector(1, 2, 3)       # 3-element column vector
w = Vector(4, 5, 6)

# Arithmetic
u = v + w                 # element-wise addition
s = 2 * v                 # scalar multiplication
dot = v @ w               # dot product

# Matrix construction (each argument is a row)
A = Matrix(Vector(1, 2), Vector(3, 4))
b = Vector(5, 6)

# Matrix-vector product
x = A @ b

# Matrix transpose and inverse
At = A.transpose()
Ainv = A.inverse()
```

### Solving linear systems

```python
from pysolverkit import LinearSystem, Vector, Matrix, LinearSolverMethod

A = Matrix(Vector(2, 1), Vector(5, 7))
b = Vector(11, 13)
system = LinearSystem(A, b)

# Direct methods
x = system.solve(LinearSolverMethod.GAUSS_ELIMINATION)

# Iterative methods
x0 = Vector(0, 0)
x = system.solve(LinearSolverMethod.GAUSS_JACOBI, tol=1e-5, initial_approximation=x0)
x = system.solve(LinearSolverMethod.GAUSS_SEIDEL, tol=1e-5, initial_approximation=x0)
```

Available methods: `LinearSolverMethod.GAUSS_ELIMINATION`, `GAUSS_JACOBI`, `GAUSS_SEIDEL`.

---

## Ordinary differential equations

### First-order IVP

Solves y′(x) = f(x, y) on [a, b] with y(a) = y₀.

```python
from pysolverkit import BivariateFunction, FirstOrderLinearODE, ODEMethod

rhs = BivariateFunction(lambda x, y: y)  # y' = y  →  solution: e^x
ode = FirstOrderLinearODE(rhs, a=0, b=1, y0=1)

# Euler's method
sol = ode.solve(h=0.1, method=ODEMethod.EULER)

# Classical 4th-order Runge-Kutta
sol = ode.solve(h=0.1, method=ODEMethod.RUNGE_KUTTA, n=4)

# Taylor series method (order n)
sol = ode.solve(h=0.1, method=ODEMethod.TAYLOR, n=2)

# Trapezoidal (implicit)
sol = ode.solve(h=0.1, method=ODEMethod.TRAPEZOIDAL)

# Adams-Bashforth (k-step explicit multistep)
sol = ode.solve(h=0.1, method=ODEMethod.ADAMS_BASHFORTH, step=4, points=[...])

# Adams-Moulton (k-step implicit multistep)
sol = ode.solve(h=0.1, method=ODEMethod.ADAMS_MOULTON, step=3, points=[...])

# Adams predictor-corrector
sol = ode.solve(h=0.1, method=ODEMethod.PREDICTOR_CORRECTOR)
```

`sol` is a `Polynomial` representing the solution; evaluate it at any point with `sol(x)`.

Available methods: `ODEMethod.EULER`, `RUNGE_KUTTA`, `TAYLOR`, `TRAPEZOIDAL`, `ADAMS_BASHFORTH`, `ADAMS_MOULTON`, `PREDICTOR_CORRECTOR`.

### Second-order BVP (linear)

Solves y″ = p(x)y′ + q(x)y + r(x) on [a, b] with boundary conditions y(a) = α, y(b) = β.

```python
from pysolverkit import (
    Function, SecondOrderLinearODE_BVP, BVPMethod
)

p = Function(lambda x: 0)
q = Function(lambda x: -4)
r = Function(lambda x: 0)

bvp = SecondOrderLinearODE_BVP(p, q, r, a=0, b=1, alpha=0, beta=2)

sol = bvp.solve(h=0.1, method=BVPMethod.SHOOTING)
sol = bvp.solve(h=0.1, method=BVPMethod.FINITE_DIFFERENCE)
```

Available methods: `BVPMethod.SHOOTING`, `FINITE_DIFFERENCE`.

### Second-order BVP (nonlinear)

Solves y″ = f(x, y, y′) with boundary conditions y(a) = α, y(b) = β.

```python
from pysolverkit import BivariateFunction, SecondOrderODE_BVP, NonlinearBVPMethod

f = BivariateFunction(lambda x, yz: yz[1]**2 - yz[0] + 1)
bvp = SecondOrderODE_BVP(f, a=0, b=1, alpha=0, beta=2)

sol = bvp.solve(h=0.1, method=NonlinearBVPMethod.SHOOTING_NEWTON, tol=1e-5)
sol = bvp.solve(h=0.1, method=NonlinearBVPMethod.FINITE_DIFFERENCE)
```

Available methods: `NonlinearBVPMethod.SHOOTING_NEWTON`, `FINITE_DIFFERENCE`.

### Second-order IVP

Reduces second-order IVPs to a first-order system.

```python
from pysolverkit import MultiVariableFunction, SecondOrderODE_IVP, ODEMethod

# y'' = -y  (simple harmonic oscillator)
f = MultiVariableFunction(lambda x, y, yp: -y)
ivp = SecondOrderODE_IVP(f, a=0, b=math.pi, y0=0, yp0=1)
sol = ivp.solve(h=0.1, method=ODEMethod.RUNGE_KUTTA, n=4)
```

---

## Utilities

`Util` provides standalone helper functions.

```python
from pysolverkit import Util

# Factorial
Util.factorial(10)

# Combinations
Util.nCr(10, 3)
```

---

## Method reference

### Root finding — `RootFindingMethod`

| Member | Algorithm | Required args |
|---|---|---|
| `BISECTION` | Bisection | `a`, `b` |
| `NEWTON` | Newton-Raphson | `p0` (derivative set via `differentiate`) |
| `MODIFIED_NEWTON` | Modified Newton | `p0` (derivative set via `differentiate`) |
| `SECANT` | Secant | `p0`, `p1` |
| `REGULA_FALSI` | Regula Falsi | `p0`, `p1` |

### Integration — `IntegrationMethod`

| Member | Rule |
|---|---|
| `RECTANGULAR` | Left-endpoint rectangle |
| `MIDPOINT` | Midpoint rule |
| `TRAPEZOIDAL` | Trapezoidal rule |
| `SIMPSON` | Simpson's 1/3 rule |
| `GAUSS` | Gauss-Legendre quadrature |

### Differentiation — `DifferentiationMethod`

| Member | Scheme |
|---|---|
| `FORWARD` | Forward difference |
| `BACKWARD` | Backward difference |
| `CENTRAL` | Central difference |

### ODE (first-order IVP) — `ODEMethod`

| Member | Method |
|---|---|
| `EULER` | Euler's method |
| `RUNGE_KUTTA` | Runge-Kutta (order n) |
| `TAYLOR` | Taylor series (order n) |
| `TRAPEZOIDAL` | Implicit trapezoidal |
| `ADAMS_BASHFORTH` | Adams-Bashforth (k-step) |
| `ADAMS_MOULTON` | Adams-Moulton (k-step) |
| `PREDICTOR_CORRECTOR` | Adams predictor-corrector |

### BVP (linear) — `BVPMethod`

| Member | Method |
|---|---|
| `SHOOTING` | Linear shooting |
| `FINITE_DIFFERENCE` | Finite difference |

### BVP (nonlinear) — `NonlinearBVPMethod`

| Member | Method |
|---|---|
| `SHOOTING_NEWTON` | Nonlinear shooting with Newton iteration |
| `FINITE_DIFFERENCE` | Finite difference |

### Linear systems — `LinearSolverMethod`

| Member | Method |
|---|---|
| `GAUSS_ELIMINATION` | Gaussian elimination |
| `GAUSS_JACOBI` | Gauss-Jacobi iteration |
| `GAUSS_SEIDEL` | Gauss-Seidel iteration |

### Interpolation — `InterpolationMethod` / `InterpolationForm`

| `InterpolationMethod` | Algorithm |
|---|---|
| `LAGRANGE` | Lagrange interpolation |
| `NEWTON` | Newton divided-difference |

| `InterpolationForm` | Form |
|---|---|
| `STANDARD` | Divided-difference table |
| `FORWARD_DIFF` | Newton forward-difference |
| `BACKWARD_DIFF` | Newton backward-difference |
