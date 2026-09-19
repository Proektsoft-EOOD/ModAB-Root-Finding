from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
import numpy as np

# IEEE-strict optimization flags per compiler. No fast-math: the solver relies on
# NaN/inf semantics, and reassociation changes the iterates. No -march=native: the
# binary must run on other CPUs. No FMA contraction, so results match other builds.
COMPILE_ARGS = {
    'msvc': ['/O2', '/fp:precise'],
    'unix': ['-O3', '-ffp-contract=off'],     # GCC and Clang
    'mingw32': ['-O3', '-ffp-contract=off'],
}


class BuildExt(build_ext):
    def build_extensions(self):
        args = COMPILE_ARGS.get(self.compiler.compiler_type, [])
        for ext in self.extensions:
            ext.extra_compile_args = args
        super().build_extensions()


extensions = [
    Extension(
        "cymodab",
        ["cymodab.pyx"],
    )
]

setup(
    name="cymodab",
    cmdclass={'build_ext': BuildExt},
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
