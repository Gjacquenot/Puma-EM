CC:= g++
MPICC:= mpiCC
# safe flags
CFLAGS:= -c -O2 -fPIC -march=native
# potentially faster flags
#CFLAGS:= -c -O2 -fPIC -ffast-math -mfpmath=both -march=native
# debug flags
#CFLAGS:= -c -g -DBZ_DEBUG -Wall -fPIC
F77:= gfortran
G2C:= gfortran
F_FLAGS:= -c -O3 -fPIC -pthread -ffast-math -mfpmath=both

# we can use vendor-supplied lapack. If empty, puma-em lapack is used
#LIBLAPACK:= -L/usr/lib/ -llapack

#INCLUDE_PATH= -I/path/to/include
#LIB_SEARCH_PATH= -L/path/to/lib
