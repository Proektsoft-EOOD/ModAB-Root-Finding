"""
Build script for the optional native extension pymodab._modab.

The extension targets the CPython Limited API (abi3), so a single binary per
platform works on all Python versions >= 3.8. If it cannot be compiled, the
package still installs and falls back to the ctypes wrapper around the shared
libraries shipped in src/pymodab.

Set PYMODAB_PURE=1 to build a pure-Python (ctypes only) py3-none-any wheel.
"""

import os
import sys
import sysconfig

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.bdist_wheel import bdist_wheel

LIMITED_API = 0x03080000
PURE = os.environ.get("PYMODAB_PURE", "") not in ("", "0")

# Prebuilt libraries for the ctypes fallback, shipped for every platform
FALLBACK_LIBS = {
    "ModAB.dll", "ModAB.lib", "ModAB.pdb",
    "libModAB.so", "libModAB_x64.dylib", "libModAB_arm64.dylib",
}


def fallback_libs_for(plat_name):
    """The fallback libraries that belong in a wheel for plat_name."""
    plat = plat_name.lower()
    if plat.startswith("win"):
        return {"ModAB.dll", "ModAB.lib"}
    if "macos" in plat:
        return {"libModAB_arm64.dylib"} if "arm64" in plat else {"libModAB_x64.dylib"}
    if "linux" in plat:
        return {"libModAB.so"}
    return set()


class BuildExt(build_ext):
    def run(self):
        super().run()
        # Once the extension is built, the fallback libraries of the other
        # platforms are dead weight. They also break wheel repair tools, which
        # reject a foreign architecture inside the wheel.
        if self.extensions and all(os.path.exists(p) for p in self.get_outputs()):
            keep = fallback_libs_for(self.plat_name)
            for lib in sorted(FALLBACK_LIBS - keep):
                path = os.path.join(self.build_lib, "pymodab", lib)
                if os.path.exists(path):
                    os.remove(path)

    def build_extensions(self):
        if self.compiler.compiler_type != "msvc":
            for ext in self.extensions:
                ext.extra_compile_args += ["-O2"]
        super().build_extensions()

    def get_libraries(self, ext):
        # With MinGW on Windows, link to python3.dll (stable ABI) instead of
        # python3XY.dll, so the abi3 binary works on every Python version.
        libs = super().get_libraries(ext)
        if sys.platform == "win32" and self.compiler.compiler_type != "msvc":
            libs = ["python3" if lib.startswith("python3") else lib for lib in libs]
        return libs


class BdistWheel(bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        if ext_modules:
            self.py_limited_api = "cp38"


ext_modules = []
# The Limited API is only available on CPython builds with the GIL
if (
    not PURE
    and sys.implementation.name == "cpython"
    and sysconfig.get_config_var("Py_GIL_DISABLED") != 1
):
    ext_modules.append(
        Extension(
            "pymodab._modab",
            ["src/pymodab/_modab.c"],
            define_macros=[("Py_LIMITED_API", hex(LIMITED_API))],
            py_limited_api=True,
            optional=True,
        )
    )

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt, "bdist_wheel": BdistWheel},
)
