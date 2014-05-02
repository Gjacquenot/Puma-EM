CC:= g++
MPICC:= mpiCC
# safe flags
CFLAGS:= -c -O2 -fPIC -ffast-math -mfpmath=both
# potentially faster flags
#CFLAGS:= -c -O3 -fPIC -ffast-math -mfpmath=both
# debug flags
#CFLAGS:= -c -g -DBZ_DEBUG -Wall -fPIC
F77:= gfortran
G2C:= gfortran
F_FLAGS:= -c -O3 -fPIC -pthread -ffast-math -mfpmath=both

#INCLUDE_PATH= -I/path/to/include
#LIB_SEARCH_PATH= -L/path/to/lib

# Do we use the lapack library furnished with the code?
# It is usually slower than vendor-supplied lapack, but
# it always works (i.e. no compilation or linking problems).
USE_VENDOR_LAPACK = 0
ifeq ($(USE_VENDOR_LAPACK), 1)
  # indicate the location and name of vendor-supplied LAPACK
  LIBLAPACK_PATH = /usr/lib/lapack/
  LIBLAPACK = lapack
else
  # default LAPACK from Puma-EM
  LIBLAPACK_PATH = ./lapack/
  LIBLAPACK = LOCAL_LAPACK
endif

