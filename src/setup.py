from setuptools import setup
from Cython.Build import cythonize
import numpy

# Build both Cython extensions in a single setup() call.  Keeping the build
# definition compact makes dependency changes easier to review and avoids two
# independent setup invocations during installation.
setup(
    name="MembIT",
    ext_modules=cythonize(["membit_module.pyx", "membrane.pyx"]),
    include_dirs=[numpy.get_include()],
    install_requires=[
        "numpy",
        "Cython",
        # MDAnalysis is only imported when reading XTC/TRR/DCD/NC trajectories,
        # but listing it here documents the optional trajectory-reader path and
        # lets standard Python tooling install a complete environment.
        "MDAnalysis>=2.0",
    ],
)
