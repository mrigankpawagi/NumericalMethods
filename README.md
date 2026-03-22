# Numerical Methods

[![Tests](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml/badge.svg)](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)

A compact Python library that brings classical numerical methods to life — root finding, interpolation, differentiation, integration, ODE solvers, and linear systems, all in a clean, composable API.

## Features

- **Function arithmetic** — add, subtract, multiply, and compose built-in function types (`Polynomial`, `Sin`, `Cos`, `Tan`, `Exponent`, `Log`) with plain Python operators
- **Root finding** — bisection, Newton-Raphson, secant, Regula Falsi, and modified Newton methods
- **Interpolation** — Lagrange and Newton divided-difference forms (forward/backward difference tables included)
- **Differentiation** — forward, backward, and central finite-difference schemes; set an exact derivative for Newton-type solvers
- **Integration** — rectangular, midpoint, trapezoidal, Simpson's, and Gauss-Legendre quadrature
- **ODE solvers** — Euler, Runge-Kutta (orders 1–4), Taylor series, trapezoidal, Adams-Bashforth/Moulton, and predictor-corrector methods for first-order IVPs; shooting and finite-difference methods for second-order BVPs
- **Linear algebra** — `Vector` and `Matrix` with arithmetic operations, and `LinearSystem` with Gaussian elimination and Gauss-Jacobi/Seidel iterative solvers

## Quick start

```python
import math
from numericalmethods import Polynomial, Sin, RootFindingMethod, IntegrationMethod

# Root finding
f = Polynomial(-6, 14, -7, 1)   # x³ − 7x² + 14x − 6
root = f.root(RootFindingMethod.BISECTION, a=0, b=1, TOLERANCE=1e-6)

# Integration
g = Sin(Polynomial(0, 1))       # sin(x)
area = g.integrate(0, math.pi, method=IntegrationMethod.SIMPSON, n=100)
```

## Documentation

Full API reference, detailed usage examples, and a method selector reference are in **[`docs/README.md`](docs/README.md)**.

## Running tests

```bash
python3 -m unittest discover -s tests -v
```

## Contributing

Contributions are welcome! To get started:

1. Fork the repository and create a feature branch.
2. Add or update tests in `tests/tests.py` for any new behaviour.
3. Ensure all tests pass: `python3 -m unittest discover -s tests -v`
4. Open a pull request with a clear description of the change.
