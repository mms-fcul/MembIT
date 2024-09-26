from distutils.core import setup
from Cython.Build import cythonize
from Cython.Compiler.Options import annotate
import numpy

#annotate = True

setup(
    ext_modules = cythonize("membit_module.pyx"),
    include_dirs=[numpy.get_include()]
)
setup(
    ext_modules = cythonize("membrane.pyx"),
    include_dirs=[numpy.get_include()]
)
