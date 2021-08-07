# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#         Make MLD - Makefile
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
default: llf.so

%.so::  %.F90
	f2py -c -m $* $<

default: brent.so

%.so::  %.F90
	f2py -c -m $* $<

default: normalize0mean.so

%.so::  %.F90
	f2py -c -m $* $<

default: pdist.so

%.so::  %.F90
	f2py -c -m $* $<