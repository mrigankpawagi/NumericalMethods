# Numerical Methods

[![Tests](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml/badge.svg)](https://github.com/mrigankpawagi/NumericalMethods/actions/workflows/tests.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)

A compact Python library for classical numerical methods: root finding, interpolation, differentiation/integration, ODE solvers, and linear systems.

## Quick start

```python
from numericalmethods import Polynomial, RootFindingMethod

f = Polynomial(-6, 14, -7, 1)
root = f.root(RootFindingMethod.BISECTION, a=0, b=1, TOLERANCE=1e-6)
```

Method selectors are typed enums (not strings), exposed from the top-level package.

## Documentation

For full API usage, package layout, enum reference, and examples, see:

- [`docs/full_documentation.md`](docs/full_documentation.md)

## Running tests

```bash
python3 -m unittest discover -s tests -v
```
