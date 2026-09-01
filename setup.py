"""Build native extensions into the wheel.

Without this, `pip wheel .` produced bwtandem-0.9.0-py3-none-any.whl: a pure-
Python wheel with no _accelerators, so installed users silently ran the slow
fallbacks while CI's manual in-place build masked the defect. The ctypes
libraries under bwtandem/c_extensions are built post-install via bwtandem/c_extensions/
build.py (they are dlopen'd by path, not imported), so only the Cython
extension belongs here.

"""
import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

setup(
    ext_modules=cythonize(
        [Extension("bwtandem._accelerators", ["bwtandem/_accelerators.pyx"],
                   include_dirs=[np.get_include()])],
        compiler_directives={"language_level": "3"},
    ),
)
