include ../../makefile.inc

LIBLAPACK_PATH = ./lapack/
LIBAMOS_ZBESH_PATH = ./amos/zbesh/
LIBITERATIVE_PATH = ./iterative/
WORKING_DIR_PATH =  $(PWD)

OBJECTS_LIBMOM = triangle_int_FS.o V_E_V_H_dipole.o V_E_V_H_plane.o integr_1D_X_W.o GL.o Z_EJ_Z_HJ_FS_triangles_arrays.o GK_triangle.o mesh.o readWriteBlitzArrayFromFile.o
OBJECTS_LIBMLFMA = triangle_int_FS.o V_E_V_H_dipole.o V_E_V_H_plane.o integr_1D_X_W.o GL.o GK_triangle.o mesh.o readWriteBlitzArrayFromFile.o interpolation.o cube.o level.o octtree.o

all: libs mpi_mlfma testMPI

testMPI: testMPI.o
	$(MPICC) testMPI.o -lblitz -lm -o testMPI

communicateMeshArrays: communicateMeshArrays.o readWriteBlitzArrayFromFile.o
	$(MPICC) communicateMeshArrays.o readWriteBlitzArrayFromFile.o -lblitz -lm -o communicateMeshArrays

communicateZnearBlocks: communicateZnearBlocks.o readWriteBlitzArrayFromFile.o
	$(MPICC) communicateZnearBlocks.o readWriteBlitzArrayFromFile.o -lblitz -lm -o communicateZnearBlocks

mpi_mlfma: mpi_mlfma.o $(OBJECTS_LIBMLFMA)
#	$(MPICC) mpi_mlfma.o -L$(WORKING_DIR_PATH) -lMoM -L$(WORKING_DIR_PATH) -lMLFMA -L$(WORKING_DIR_PATH)/amos/zbesh -lAMOS -l$(G2C) -lblitz -lm -o mpi_mlfma
	$(MPICC) mpi_mlfma.o -L$(WORKING_DIR_PATH) -lMLFMA -L$(WORKING_DIR_PATH)/amos/zbesh -lAMOS -L$(WORKING_DIR_PATH)/lapack -lZGELS -l$(G2C) -lblitz -lm -o mpi_mlfma

libs: libMoM libMLFMA libAMOS libITERATIVE libLAPACK

libMoM: $(OBJECTS_LIBMOM)
	ar -r libMoM.a $(OBJECTS_LIBMOM)

libMLFMA: $(OBJECTS_LIBMLFMA)
	ar -r libMLFMA.a $(OBJECTS_LIBMLFMA)

libAMOS:
	cd $(LIBAMOS_ZBESH_PATH); make libAMOS

libITERATIVE:
	cd $(LIBITERATIVE_PATH); make libITERATIVE

libLAPACK:
	cd $(LIBLAPACK_PATH); make

testMPI.o: testMPI.cpp
	$(MPICC) $(CFLAGS) testMPI.cpp
communicateMeshArrays.o: communicateMeshArrays.cpp readWriteBlitzArrayFromFile.h
	$(MPICC) $(CFLAGS) communicateMeshArrays.cpp
communicateZnearBlocks.o: communicateZnearBlocks.cpp readWriteBlitzArrayFromFile.h
	$(MPICC) $(CFLAGS) communicateZnearBlocks.cpp
mpi_mlfma.o: mpi_mlfma.cpp  Z_sparse_MLFMA.h mesh.h octtree.h octtree.cpp level.h level.cpp ./iterative/iterative.h V_E_V_H.h alpha_computation.h V_E_V_H_dipole.cpp V_E_V_H_plane.cpp EMConstants.h
	$(MPICC) $(CFLAGS) mpi_mlfma.cpp
octtree.o: octtree.cpp interpolation.h level.h level.cpp FMM.h mesh.h octtree.h integr_1D_X_W.h alpha_computation.h readWriteBlitzArrayFromFile.h
	$(MPICC) $(CFLAGS) octtree.cpp
level.o: level.cpp interpolation.h level.h FMM.h mesh.h cube.h cube.cpp integr_1D_X_W.h readWriteBlitzArrayFromFile.h alpha_computation.h
	$(MPICC) $(CFLAGS) level.cpp
cube.o: cube.cpp FMM.h mesh.h cube.h integr_1D_X_W.h readWriteBlitzArrayFromFile.h
	$(MPICC) $(CFLAGS) cube.cpp
FMM.o: FMM.cpp EMConstants.h FMM.h integr_1D_X_W.h GK_triangle.h $(LIBAMOS_ZBESH_PATH)zbesh_interface.h
	$(CC) $(CFLAGS) FMM.cpp
integr_1D_X_W.o: integr_1D_X_W.cpp GL.h 
	$(CC) $(CFLAGS) integr_1D_X_W.cpp
interpolation.o: interpolation.cpp interpolation.h
	$(CC) $(CFLAGS) interpolation.cpp
mesh.o: mesh.cpp mesh.h
	$(CC) $(CFLAGS) mesh.cpp
readWriteBlitzArrayFromFile.o: readWriteBlitzArrayFromFile.cpp readWriteBlitzArrayFromFile.h
	$(CC) $(CFLAGS) readWriteBlitzArrayFromFile.cpp
Z_EJ_Z_HJ_FS_triangles_arrays.o: Z_EJ_Z_HJ_FS_triangles_arrays.cpp EMConstants.h triangle_int.h mesh.h
	$(CC) $(CFLAGS) Z_EJ_Z_HJ_FS_triangles_arrays.cpp
V_E_V_H_dipole.o: V_E_V_H_dipole.cpp GK_triangle.h triangle_int.h EMConstants.h V_E_V_H.h
	$(CC) $(CFLAGS) V_E_V_H_dipole.cpp
V_E_V_H_plane.o: V_E_V_H_plane.cpp GK_triangle.h triangle_int.h EMConstants.h V_E_V_H.h
	$(CC) $(CFLAGS) V_E_V_H_plane.cpp
triangle_int_FS.o: triangle_int_FS.cpp GK_triangle.h GL.h
	$(CC) $(CFLAGS) triangle_int_FS.cpp
GK_triangle.o: GK_triangle.cpp GK_triangle.h
	$(CC) $(CFLAGS) GK_triangle.cpp
GL.o: GL.cpp
	$(CC) $(CFLAGS) GL.cpp



clean: 
	rm -f *.o *.a *.txt *~ *.pyc octtree mpi_mlfma testMPI communicateZnearBlocks communicateMeshArrays; 
	cd $(LIBAMOS_ZBESH_PATH); make clean; 
	cd $(LIBLAPACK_PATH); make clean; 
	cd $(LIBITERATIVE_PATH); make clean;

