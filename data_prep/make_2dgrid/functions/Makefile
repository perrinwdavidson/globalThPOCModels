# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#     Make 2D Model Grid - Makefile
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
default: interp2d.so

%.so::  %.F90
	f2py -c -m $* $<