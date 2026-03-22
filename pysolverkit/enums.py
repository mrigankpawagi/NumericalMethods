from enum import Enum


class DifferentiationMethod(Enum):
    """Numerical differentiation method."""

    FORWARD = "forward"
    BACKWARD = "backward"
    CENTRAL = "central"


class IntegrationMethod(Enum):
    """Numerical integration (quadrature) method."""

    RECTANGULAR = "rectangular"
    MIDPOINT = "midpoint"
    TRAPEZOIDAL = "trapezoidal"
    SIMPSON = "simpson"
    GAUSS = "gauss"


class RootFindingMethod(Enum):
    """Root-finding method for scalar equations."""

    BISECTION = "bisection"
    NEWTON = "newton"
    MODIFIED_NEWTON = "modified_newton"
    SECANT = "secant"
    REGULA_FALSI = "regula_falsi"


class InterpolationMethod(Enum):
    """Polynomial interpolation method."""

    LAGRANGE = "lagrange"
    NEWTON = "newton"


class InterpolationForm(Enum):
    """Form of Newton interpolating polynomial."""

    STANDARD = "standard"
    FORWARD_DIFF = "forward_diff"
    BACKWARD_DIFF = "backward_diff"


class ODEMethod(Enum):
    """Numerical solver for first-order ODE initial value problems."""

    EULER = "euler"
    RUNGE_KUTTA = "runge-kutta"
    TAYLOR = "taylor"
    TRAPEZOIDAL = "trapezoidal"
    ADAMS_BASHFORTH = "adam-bashforth"
    ADAMS_MOULTON = "adam-moulton"
    PREDICTOR_CORRECTOR = "predictor-corrector"


class BVPMethod(Enum):
    """Solver for second-order linear boundary value problems."""

    SHOOTING = "shooting"
    FINITE_DIFFERENCE = "finite_difference"


class NonlinearBVPMethod(Enum):
    """Solver for second-order nonlinear boundary value problems."""

    SHOOTING_NEWTON = "shooting_newton"
    FINITE_DIFFERENCE = "finite_difference"


class LinearSolverMethod(Enum):
    """Solver for systems of linear equations."""

    GAUSS_ELIMINATION = "gauss_elimination"
    GAUSS_JACOBI = "gauss_jacobi"
    GAUSS_SEIDEL = "gauss_seidel"


__all__ = [
    "DifferentiationMethod",
    "IntegrationMethod",
    "RootFindingMethod",
    "InterpolationMethod",
    "InterpolationForm",
    "ODEMethod",
    "BVPMethod",
    "NonlinearBVPMethod",
    "LinearSolverMethod",
]
