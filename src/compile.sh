for i in 6 7 8 9 10
do
    # python3.$i -m pip install cython numpy
    python3.$i setup.py build_ext --inplace
done
