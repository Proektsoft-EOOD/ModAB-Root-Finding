from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "cymodab",
        ["cymodab.pyx"],
        extra_compile_args=['-O3', '-march=native', '-ffast-math'],
        extra_link_args=[],
    )
]

setup(
    name="cymodab",
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            'language_level': "3",
            'boundscheck': False,
            'wraparound': False,
            'cdivision': True,
            'initializedcheck': False,
        },
        annotate=True,  # Creates HTML file showing Python/C interaction
    ),
)
