#!/bin/bash

cd generalized_zeta_leskovec_code/main
make clean
rm gen_zeta.so
rm gen_zeta.exe
make
