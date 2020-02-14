from distutils.core import setup
from Cython.Build import cythonize
from Cython.Compiler.Options import annotate

#annotate = True

setup(
    ext_modules = cythonize("membit_module.pyx")
)
setup(
    ext_modules = cythonize("membrane.pyx")
)
