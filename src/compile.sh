#!/bin/bash -e

PYTHON_BIN="${PYTHON_BIN:-python}"

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    PYTHON_BIN="python3"
fi

# Use the active Python environment instead of a hard-coded interpreter.  This
# makes it safe to run inside a virtualenv/conda env dedicated to XTC testing.
"${PYTHON_BIN}" -m pip install --upgrade cython numpy MDAnalysis
"${PYTHON_BIN}" setup.py build_ext --inplace
