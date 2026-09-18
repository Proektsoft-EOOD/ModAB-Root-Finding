"""
Python wrapper for the ModAB C library.

Uses the native extension module pymodab._modab when it is available, and falls
back to calling the bundled shared library through ctypes otherwise.
"""

import ctypes
import os
import platform
import sys
from ctypes import c_double, c_int, CFUNCTYPE
from typing import Callable, Union

# Define the callback function type for the C library
FUNC_TYPE = CFUNCTYPE(c_double, c_double)

def _get_library_path() -> str:
    """Get the path to the shared library based on the platform."""
    package_dir = os.path.dirname(os.path.abspath(__file__))

    system = platform.system()
    if system == "Windows":
        lib_name = "ModAB.dll"
    elif system == "Darwin":
        machine = platform.machine()
        if machine == "arm64":
            lib_name = "libModAB_arm64.dylib"
        else:
            lib_name = "libModAB_x64.dylib"
    else:  # Linux and others
        lib_name = "libModAB.so"

    lib_path = os.path.join(package_dir, lib_name)

    if not os.path.exists(lib_path):
        raise FileNotFoundError(
            f"Could not find {lib_name} in {package_dir}. "
            f"Please ensure the library is built for your platform."
        )

    return lib_path

def _load_library():
    """Load the ModAB shared library."""
    lib_path = _get_library_path()
    lib = ctypes.CDLL(lib_path)

    # Configure modAB_find_root
    lib.modAB_find_root.argtypes = [FUNC_TYPE, c_double, c_double, c_double, c_double, c_int]
    lib.modAB_find_root.restype = c_double

    # Configure get_evaluation_count
    lib.get_evaluation_count.argtypes = []
    lib.get_evaluation_count.restype = c_int

    return lib


def _ctypes_find_root(
    f: Union[Callable[[float], float], int],
    x1: float,
    x2: float,
    atol: float = 1e-14,
    rtol: float = 1e-14,
    max_iter: int = 200
) -> float:
    """
    Find the root of f(x) = 0 within the interval [x1, x2].

    Uses an improved version of the Modified Anderson-Bjork method
    (Ganchovski, Traykov) which combines bisection with the secant method
    for robust and fast convergence.

    Parameters
    ----------
    f : callable or int
        A continuous function of one variable, or the address (int) of a
        compiled C function ``double f(double)``, e.g. ``numba.cfunc(...).address``.
        A compiled function is called directly, without any Python overhead.
    x1 : float
        Left endpoint of the bracket interval.
    x2 : float
        Right endpoint of the bracket interval.
    atol : float, optional
        Absolute tolerance for convergence (default: 1e-14).
    rtol : float, optional
        Relative tolerance for convergence (default: 1e-14).
    max_iter : int, optional
        Maximum number of iterations (default: 200).

    Returns
    -------
    float
        The root of f(x) = 0 within [x1, x2].
        Returns NaN if no root is found or if f(x1) and f(x2) have the same sign.

    Notes
    -----
    The function f must be continuous on [x1, x2] and f(x1) * f(x2) < 0
    (i.e., the function must have opposite signs at the endpoints).
    Exceptions raised by f are propagated to the caller (native extension only).

    Examples
    --------
    >>> import math
    >>> from pymodab import find_root
    >>> # Find the root of cos(x) - x in [0, 1]
    >>> root = find_root(lambda x: math.cos(x) - x, 0, 1)
    >>> print(f"{root:.10f}")
    0.7390851332
    """
    c_func = FUNC_TYPE(f)
    return _lib.modAB_find_root(c_func, x1, x2, atol, rtol, max_iter)


def _ctypes_get_evaluation_count() -> int:
    """
    Get the number of function evaluations from the last root-finding call.

    Returns
    -------
    int
        The number of times the function was evaluated during the last
        call to find_root.

    Examples
    --------
    >>> from pymodab import find_root, get_evaluation_count
    >>> root = find_root(lambda x: x**2 - 2, 1, 2)
    >>> print(f"Evaluations: {get_evaluation_count()}")
    Evaluations: 11
    """
    return _lib.get_evaluation_count()


try:
    from . import _modab
except ImportError:
    _modab = None

NATIVE = _modab is not None
"""True if the native extension is used, False if the ctypes fallback is."""

if NATIVE:
    # Bind the C functions directly, so no Python code runs per call
    find_root = _modab.find_root
    get_evaluation_count = _modab.get_evaluation_count
else:
    _lib = _load_library()
    find_root = _ctypes_find_root
    get_evaluation_count = _ctypes_get_evaluation_count
