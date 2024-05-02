# Numerical Methods and Functions Library

**Created during UMC 202 at the Indian Institute of Science, Bengaluru**

<details>
  <summary><b>Functions</b></summary>

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

**Other Utilities**
- Plotting (requires `matplotlib`)
</details>

<details>
  <summary><b>Methods</b></summary>

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
</details>

## Usage

The file `function.py` contains all the important classes and methods, and is the main file to be imported. The file `util.py` has a few helper methods required by the main file. Some practice problem sets have been included in the `examples` directory to demonstrate the usage of the library.

> [!NOTE]  
> Documentation below was generated automatically using [pdoc](https://pdoc3.github.io/).

# Module `function` 








## Classes



### Class `BivariateFunction` 




>     class BivariateFunction(
>         function
>     )






#### Ancestors (in MRO)

* [function.MultiVariableFunction](#function.MultiVariableFunction)







### Class `Cos` 




>     class Cos(
>         f: function.Function
>     )






#### Ancestors (in MRO)

* [function.Function](#function.Function)







### Class `Exponent` 




>     class Exponent(
>         f: function.Function,
>         base: float = 2.718281828459045
>     )






#### Ancestors (in MRO)

* [function.Function](#function.Function)







### Class `FirstOrderLinearODE` 




>     class FirstOrderLinearODE(
>         f: function.BivariateFunction,
>         a: float,
>         b: float,
>         y0: float
>     )


y'(x) = f(x, y(x))
These are initial value problems.

f is a function of x and y(x)



#### Ancestors (in MRO)

* [function.LinearODE](#function.LinearODE)
* [function.OrdinaryDifferentialEquation](#function.OrdinaryDifferentialEquation)







#### Methods



##### Method `solve` 




>     def solve(
>         self,
>         h: float = 0.1,
>         method: Literal['euler', 'runge-kutta', 'taylor', 'trapezoidal', 'adam-bashforth', 'adam-moulton', 'predictor-corrector'] = 'euler',
>         n: int = 1,
>         step: int = 2,
>         points: list[float] = []
>     )





##### Method `solve_adam_bashforth` 




>     def solve_adam_bashforth(
>         self,
>         h: float,
>         step: int,
>         points: list[float]
>     ) ‑> function.Polynomial





##### Method `solve_adam_moulton` 




>     def solve_adam_moulton(
>         self,
>         h: float,
>         step: int,
>         points: list[float]
>     ) ‑> function.Polynomial





##### Method `solve_predictor_corrector` 




>     def solve_predictor_corrector(
>         self,
>         h: float
>     ) ‑> function.Polynomial





##### Method `solve_runge_kutta` 




>     def solve_runge_kutta(
>         self,
>         h: float,
>         n: int
>     ) ‑> function.Polynomial





##### Method `solve_taylor` 




>     def solve_taylor(
>         self,
>         h: float,
>         n: int
>     ) ‑> function.Polynomial





##### Method `solve_trapezoidal` 




>     def solve_trapezoidal(
>         self,
>         h: float
>     ) ‑> function.Polynomial





### Class `Function` 




>     class Function(
>         function
>     )







#### Descendants

* [function.Cos](#function.Cos)
* [function.Exponent](#function.Exponent)
* [function.Log](#function.Log)
* [function.Polynomial](#function.Polynomial)
* [function.Sin](#function.Sin)
* [function.Tan](#function.Tan)






#### Methods



##### Method `bisection` 




>     def bisection(
>         self,
>         a: float,
>         b: float,
>         TOLERANCE=1e-10,
>         N=100,
>         early_stop: int = None
>     )





##### Method `differentiate` 




>     def differentiate(
>         self,
>         func=None,
>         h=1e-05,
>         method: Literal['forward', 'backward', 'central'] = 'forward'
>     )


Sets or returns the derivative of the function.
If func is None, returns the derivative.
If func is a Function, sets the derivative to func.
If func is a lambda, sets the derivative to a Function with the lambda.


##### Method `differentiate_central` 




>     def differentiate_central(
>         self,
>         h
>     )





##### Method `differentiate_forward` 




>     def differentiate_forward(
>         self,
>         h
>     )





##### Method `fixed_point` 




>     def fixed_point(
>         self,
>         p0: float,
>         TOLERANCE=1e-10,
>         N=100
>     )





##### Method `integral` 




>     def integral(
>         self,
>         func=None,
>         h=1e-05
>     )


Sets or returns the integral of the function.
If func is None, returns the integral.
If func is a Function, sets the integral to func.
If func is a lambda, sets the integral to a Function with the lambda.


##### Method `integrate` 




>     def integrate(
>         self,
>         a: float,
>         b: float,
>         method: Literal['rectangular', 'midpoint', 'trapezoidal', 'simpson', 'gauss'] = None,
>         n: int = None
>     )


Definite integral of the function from a to b.


##### Method `integrate_gauss` 




>     def integrate_gauss(
>         self,
>         a: float,
>         b: float,
>         n: int = None
>     )





##### Method `integrate_midpoint` 




>     def integrate_midpoint(
>         self,
>         a: float,
>         b: float,
>         n: int = None
>     )





##### Method `integrate_rectangular` 




>     def integrate_rectangular(
>         self,
>         a: float,
>         b: float,
>         n: int = None
>     )





##### Method `integrate_simpson` 




>     def integrate_simpson(
>         self,
>         a: float,
>         b: float,
>         n: int = None
>     )





##### Method `integrate_trapezoidal` 




>     def integrate_trapezoidal(
>         self,
>         a: float,
>         b: float,
>         n: int = None
>     )





##### Method `modified_newton` 




>     def modified_newton(
>         self,
>         p0: float,
>         TOLERANCE=1e-10,
>         N=100,
>         early_stop: int = None
>     )





##### Method `multi_differentiate` 




>     def multi_differentiate(
>         self,
>         n: int,
>         h=1e-05,
>         method: Literal['forward', 'backward', 'central'] = 'forward'
>     )


Returns the nth derivative of the function.


##### Method `newton` 




>     def newton(
>         self,
>         p0: float,
>         TOLERANCE=1e-10,
>         N=100,
>         early_stop: int = None
>     )





##### Method `plot` 




>     def plot(
>         self,
>         min: float,
>         max: float,
>         N=1000,
>         file: str = '',
>         clear: bool = False
>     )





##### Method `regula_falsi` 




>     def regula_falsi(
>         self,
>         p0: float,
>         p1: float,
>         TOLERANCE=1e-10,
>         N=100,
>         early_stop: int = None
>     )





##### Method `root` 




>     def root(
>         self,
>         method: Literal['bisection', 'newton', 'secant', 'regula_falsi', 'modified_newton'],
>         a: float = None,
>         b: float = None,
>         p0: float = None,
>         p1: float = None,
>         TOLERANCE=1e-10,
>         N=100,
>         return_iterations=False,
>         early_stop: int = None
>     )





##### Method `secant` 




>     def secant(
>         self,
>         p0: float,
>         p1: float,
>         TOLERANCE=1e-10,
>         N=100,
>         early_stop: int = None
>     )





### Class `LinearODE` 




>     class LinearODE






#### Ancestors (in MRO)

* [function.OrdinaryDifferentialEquation](#function.OrdinaryDifferentialEquation)



#### Descendants

* [function.FirstOrderLinearODE](#function.FirstOrderLinearODE)
* [function.SecondOrderLinearODE_BVP](#function.SecondOrderLinearODE_BVP)






### Class `LinearSystem` 




>     class LinearSystem(
>         A: function.Matrix,
>         b: function.Vector
>     )


A system of linear equations.








#### Methods



##### Method `solve` 




>     def solve(
>         self,
>         method: Literal['gauss_elimination', 'gauss_jacobi', 'gauss_seidel'] = 'gauss_elimination',
>         TOL: float = 1e-05,
>         initial_approximation: function.Vector = None,
>         MAX_ITERATIONS: int = 100
>     )





##### Method `solve_gauss_elimination` 




>     def solve_gauss_elimination(
>         self
>     ) ‑> function.Vector





##### Method `solve_gauss_jacobi` 




>     def solve_gauss_jacobi(
>         self,
>         TOL: float,
>         initial_approximation: function.Vector,
>         MAX_ITERATIONS
>     ) ‑> function.Vector





##### Method `solve_gauss_seidel` 




>     def solve_gauss_seidel(
>         self,
>         TOL: float,
>         initial_approximation: function.Vector,
>         MAX_ITERATIONS
>     ) ‑> function.Vector





### Class `Log` 




>     class Log(
>         f: function.Function,
>         base: float = 2.718281828459045
>     )






#### Ancestors (in MRO)

* [function.Function](#function.Function)







### Class `Matrix` 




>     class Matrix(
>         *rows: list[function.Vector]
>     )











### Class `MultiVariableFunction` 




>     class MultiVariableFunction(
>         function
>     )







#### Descendants

* [function.BivariateFunction](#function.BivariateFunction)






### Class `OrdinaryDifferentialEquation` 




>     class OrdinaryDifferentialEquation







#### Descendants

* [function.LinearODE](#function.LinearODE)
* [function.SecondOrderODE_BVP](#function.SecondOrderODE_BVP)
* [function.SecondOrderODE_IVP](#function.SecondOrderODE_IVP)






### Class `Polynomial` 




>     class Polynomial(
>         *coefficients
>     )


coefficients are in the form a_0, a_1, ... a_n



#### Ancestors (in MRO)

* [function.Function](#function.Function)






#### Static methods



##### `Method interpolate` 




>     def interpolate(
>         data: tuple,
>         method: Literal['lagrange', 'newton'] = 'newton',
>         f: function.Function = None,
>         form: Literal['standard', 'backward_diff', 'forward_diff'] = 'standard'
>     )


data is a list of (x, y) tuples.
alternative: f is a Function that returns the y values and data is a list of x values.


##### `Method interpolate_lagrange` 




>     def interpolate_lagrange(
>         data: tuple
>     )


data is a tuple of (x, y) tuples


##### `Method interpolate_newton` 




>     def interpolate_newton(
>         data: tuple
>     )


data is a tuple of (x, y) tuples


##### `Method interpolate_newton_backward_diff` 




>     def interpolate_newton_backward_diff(
>         data: tuple
>     )


data is a tuple of (x, y) tuples


##### `Method interpolate_newton_forward_diff` 




>     def interpolate_newton_forward_diff(
>         data: tuple
>     )


data is a tuple of (x, y) tuples



### Class `SecondOrderLinearODE_BVP` 




>     class SecondOrderLinearODE_BVP(
>         p: function.Function,
>         q: function.Function,
>         r: function.Function,
>         a: float,
>         b: float,
>         y0: float,
>         y1: float
>     )


y''(x) = p(x)y'(x) + q(x)y(x) + r(x)
These are boundary value problems.



#### Ancestors (in MRO)

* [function.LinearODE](#function.LinearODE)
* [function.OrdinaryDifferentialEquation](#function.OrdinaryDifferentialEquation)







#### Methods



##### Method `solve` 




>     def solve(
>         self,
>         h: float = 0.1,
>         method: Literal['shooting', 'finite_difference'] = 'shooting'
>     )





##### Method `solve_finite_difference` 




>     def solve_finite_difference(
>         self,
>         h: float
>     ) ‑> function.Polynomial





##### Method `solve_shooting` 




>     def solve_shooting(
>         self,
>         h: float
>     ) ‑> function.Polynomial





### Class `SecondOrderODE_BVP` 




>     class SecondOrderODE_BVP(
>         f: function.MultiVariableFunction,
>         a: float,
>         b: float,
>         y0: float,
>         y1: float
>     )


y''(x) = f(x, y(x), y'(x))
These are boundary value problems.

f is a function of x, y(x), and y'(x)



#### Ancestors (in MRO)

* [function.OrdinaryDifferentialEquation](#function.OrdinaryDifferentialEquation)







#### Methods



##### Method `solve` 




>     def solve(
>         self,
>         h: float = 0.1,
>         method: Literal['shooting_newton'] = 'shooting_newton',
>         M: int = 100,
>         TOL: float = 1e-05,
>         initial_approximation=None
>     )





##### Method `solve_shooting_newton` 




>     def solve_shooting_newton(
>         self,
>         h: float,
>         M,
>         TOL,
>         initial_approximation
>     ) ‑> function.Polynomial





### Class `SecondOrderODE_IVP` 




>     class SecondOrderODE_IVP(
>         f: function.MultiVariableFunction,
>         a: float,
>         b: float,
>         y0: float,
>         y1: float
>     )


y''(x) = f(x, y(x), y'(x))
These are initial value problems.

f is a function of x, y(x), and y'(x)



#### Ancestors (in MRO)

* [function.OrdinaryDifferentialEquation](#function.OrdinaryDifferentialEquation)







#### Methods



##### Method `solve` 




>     def solve(
>         self,
>         h: float = 0.1,
>         method: Literal['euler', 'runge-kutta', 'taylor', 'trapezoidal', 'adam-bashforth', 'adam-moulton', 'predictor-corrector'] = 'euler',
>         n: int = 1,
>         step: int = 2,
>         points: list[float] = []
>     )





### Class `Sin` 




>     class Sin(
>         f: function.Function
>     )






#### Ancestors (in MRO)

* [function.Function](#function.Function)







### Class `Tan` 




>     class Tan(
>         f: function.Function
>     )






#### Ancestors (in MRO)

* [function.Function](#function.Function)







### Class `Vector` 




>     class Vector(
>         *components
>     )











