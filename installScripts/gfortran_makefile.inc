CC:= g++
MPICC:= mpiCC
CFLAGS:= -c -O2 -fPIC -ffast-math
#CFLAGS:= -c -O2 -fPIC -ffast-math -ftree-vectorize -msse2
#CFLAGS:= -c -g -DBZ_DEBUG -fPIC
F77:= gfortran
G2C:= gfortran
F_FLAGS:= -c -O3 -funroll-loops -fPIC -pthread

#INCLUDE_PATH= -I/path/to/include
#LIB_SEARCH_PATH= -L/path/to/lib

