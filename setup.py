"""Build native extensions into the wheel.

Without this, `pip wheel .` produced bwtandem-0.9.0-py3-none-any.whl: a pure-
Python wheel with no _accelerators, so installed users silently ran the slow
fallbacks while CI's manual in-place build masked the defect. The ctypes
libraries under src/c_extensions are built post-install via src/c_extensions/
build.py (they are dlopen'd by path, not imported), so only the Cython
extension belongs here.

NOTE (pre-release): the import package is still `src`, which is collision-
prone for an installed distribution; renaming it to `bwtandem` is queued for
the release commit (#14) and will change the entry point to
`bwtandem.main:main`.
"""
import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

setup(
    ext_modules=cythonize(
        [Extension("src._accelerators", ["src/_accelerators.pyx"],
                   include_dirs=[np.get_include()])],
        compiler_directives={"language_level": "3"},
    ),
)
