# Full Documentation

## Package layout

```
numericalmethods/
├── enums.py              # Typed enums for method selectors
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

## Enums

- `DifferentiationMethod`: `FORWARD`, `BACKWARD`, `CENTRAL`
- `IntegrationMethod`: `RECTANGULAR`, `MIDPOINT`, `TRAPEZOIDAL`, `SIMPSON`, `GAUSS`
- `RootFindingMethod`: `BISECTION`, `NEWTON`, `MODIFIED_NEWTON`, `SECANT`, `REGULA_FALSI`
- `InterpolationMethod`: `LAGRANGE`, `NEWTON`
- `InterpolationForm`: `STANDARD`, `FORWARD_DIFF`, `BACKWARD_DIFF`
- `ODEMethod`: `EULER`, `RUNGE_KUTTA`, `TAYLOR`, `TRAPEZOIDAL`, `ADAMS_BASHFORTH`, `ADAMS_MOULTON`, `PREDICTOR_CORRECTOR`
- `BVPMethod`: `SHOOTING`, `FINITE_DIFFERENCE`
- `NonlinearBVPMethod`: `SHOOTING_NEWTON`, `FINITE_DIFFERENCE`
- `LinearSolverMethod`: `GAUSS_ELIMINATION`, `GAUSS_JACOBI`, `GAUSS_SEIDEL`

## Usage examples

```python
from numericalmethods import (
    Polynomial, Sin, BivariateFunction,
    FirstOrderLinearODE, LinearSystem, Vector, Matrix,
    RootFindingMethod, IntegrationMethod, ODEMethod, LinearSolverMethod,
)

# Root finding
f = Polynomial(-6, 14, -7, 1)
root = f.root(RootFindingMethod.BISECTION, a=0, b=1, TOLERANCE=1e-6)

# Integration
g = Sin(Polynomial(0, 1))
area = g.integrate(0, 3.14159, method=IntegrationMethod.SIMPSON, n=100)

# ODE
rhs = BivariateFunction(lambda x, y: y)
ode = FirstOrderLinearODE(rhs, a=0, b=1, y0=1)
sol = ode.solve(h=0.1, method=ODEMethod.RUNGE_KUTTA, n=4)

# Linear system
A = Matrix(Vector(2, 1), Vector(5, 7))
b = Vector(11, 13)
x = LinearSystem(A, b).solve(LinearSolverMethod.GAUSS_ELIMINATION)
```
