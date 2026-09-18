@echo off
REM Build sdist and wheels for PyPI distribution:
REM   pymodab-*-cp38-abi3-win_amd64.whl - native extension (Windows, Python 3.8+)
REM   pymodab-*-py3-none-any.whl        - pure ctypes fallback for all other platforms
REM Without MSVC, MinGW gcc can be used by setting DIST_EXTRA_CONFIG to a file containing:
REM   [build_ext]
REM   compiler=mingw32
pip install build
py -m build
set PYMODAB_PURE=1
py -m build --wheel
set PYMODAB_PURE=
echo.
echo Wheels created in dist/ directory
echo Upload to PyPI with: twine upload dist/*
