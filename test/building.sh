#! /bin/sh
gfortran -o libreentry_sub.so sub_reentry.f90 -shared -cpp -fPIC -I ../build/ -L../build -lreentry_calc
