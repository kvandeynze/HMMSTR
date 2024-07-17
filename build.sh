#!/bin/bash

chmod +x $SRC_DIR/src/c_files/test_hmm_cython.py
cp $SRC_DIR/src/c_files/test_hmm_cython.py $PREFIX/bin
$PYTHON setup.py clean
$PYTHON setup.py build_ext
$PYTHON setup.py install --single-version-externally-managed --record=record.txt