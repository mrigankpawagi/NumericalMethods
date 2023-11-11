# Numerical Methods and Functions Library

## Created during UMC 202 at the Indian Institute of Science, Bengaluru

### Functions
- Support for polynomials
- Support for function arithmetic, composition and other methods
- Support for derivatives and differentiation
    - Forward difference
    - Backward difference
    - Central difference
    - n-th order derivative
- Allows probing algorithms
- Support for integration
    - Rectangular
    - Midpoint
    - Trapezoidal 
    - Simpson's
    - Gaussian Quadrature (n = 1, 2)
- Support for multivariable functions
- Support for vectors and vector-valued functions
- Support for matrices (vector of vectors)

### Methods
- Root finding
    - Bisection
    - Newton-Raphson and Modified Newton
    - Fixed-point iteration
    - Secant
    - Regula-Falsi
- Interpolating Polynomial
    - Lagrange
    - Newton's divided difference
    - Forward difference
    - Backward difference
- Solving initial value problems on one-dimensional first order linear ODEs (and systems of such ODEs)
    - Euler's method
    - Taylor's method (for n = 1, 2)
    - Runga Kutta (for n = 1, 2, 3, 4)
    - Trapezoidal
    - Adam-Bashforth (for n = 2, 3, 4)
    - Adam-Moulton (for n = 2, 3, 4)
    - Predictor-Corrector (with initial approximation from Runga-Kutta order-4, predictor as Adam-Bashforth order-4 and corrector as Adam-Moulton order-3)
- Solving boundary value problems on one-dimensional second order linear ODEs
    - Shooting method
    - Finite difference method
- Solving initial value problems on one-dimensional second order ODEs
    - All methods from solving first order linear ODEs
- Solving boundary value problems on one-dimensional second order ODEs
    - Shooting method
        - Newton method
- Solving systems of linear equations
    - Gaussian elimination with backward substitution
    - Gauss-Jacobi method
    - Gauss-Seidel method

### Utilities
- Plotting (requires `matplotlib`)

---
Some practice problem sets have been included in the `examples` directory.
