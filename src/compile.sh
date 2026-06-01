#!/bin/bash -e

python -m pip install --upgrade cython numpy
python setup.py build_ext --inplace
