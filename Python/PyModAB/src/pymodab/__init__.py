"""
PyModAB - Python wrapper for the Modified Anderson-Bjork (modAB) algorithm .
ModAB is a fast and robust root-finding method that combines bisection with the 
Anderson-Bjork method for that provide fast convergence for well-behaved functions
while preserving wost-case optimality.

Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. 
Improvements of the Modified Anderson-Björck (modAB) Root-Finding Algorithm. 
Preprints 2026, 2026032190. https://doi.org/10.20944/preprints202603.2190.v1
"""

from .modab import find_root, get_evaluation_count, NATIVE

__version__ = "1.0.10"
__all__ = ["find_root", "get_evaluation_count"]
