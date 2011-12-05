include ../../makefile.inc

LIBLAPACK_PATH = ./lapack/
LIBAMOS_ZBESH_PATH = ./amos/zbesh/
LIBITERATIVE_PATH = ./iterative/
WORKING_DIR_PATH =  $(PWD)

OBJECTS_LIBMOM = triangle_int_FS.o V_E_V_H_dipole.o V_E_V_H_plane.o integr_1D_X_W.o GL.o Z_EJ_Z_HJ_FS_triangles_arrays.o GK_triangle.o mesh.o readWriteBlitzArrayFromFile.o
OBJECTS_LIBMLFMA = triangle_int_FS.o V_E_V_H_dipole.o V_E_V_H_plane.o integr_1D_X_W.o GL.o GK_triangle.o mesh.o readWriteBlitzArrayFromFile.o interpolation.o cube.o level.o octtree.o

all: libs mpi_mlfma testMPI

testMPI: testMPI.o
	$(MPICC) $(INCLUDE_PATH) testMPI.o $(LIB_SEARCH_PATH) -lblitz -lm -o testMPI

communicateZnearBlocks: communicateZnearBlocks.o readWriteBlitzArrayFromFile.o
	$(MPICC) $(INCLUDE_PATH) communicateZnearBlocks.o readWriteBlitzArrayFromFile.o $(LIB_SEARCH_PATH) -lblitz -lm -o communicateZnearBlocks

RWGs_renumbering: RWGs_renumbering.o readWriteBlitzArrayFromFile.o
	$(MPICC) $(INCLUDE_PATH) RWGs_renumbering.o readWriteBlitzArrayFromFile.o $(LIB_SEARCH_PATH) -lblitz -lm -o RWGs_renumbering

scatter_mesh_per_cube: scatter_mesh_per_cube.o readWriteBlitzArrayFromFile.o
	$(MPICC) $(INCLUDE_PATH) scatter_mesh_per_cube.o readWriteBlitzArrayFromFile.o $(LIB_SEARCH_PATH) -lblitz -lm -o scatter_mesh_per_cube

compute_Z_near: compute_Z_near.o $(OBJECTS_LIBMOM)
	$(MPICC) $(INCLUDE_PATH) compute_Z_near.o -L$(WORKING_DIR_PATH) -lMoM $(LIB_SEARCH_PATH) -lblitz -lm -o compute_Z_near

mpi_mlfma: mpi_mlfma.o $(OBJECTS_LIBMLFMA)
#	$(MPICC) $(INCLUDE_PATH) mpi_mlfma.o -L$(WORKING_DIR_PATH) -lMoM -L$(WORKING_DIR_PATH) -lMLFMA -L$(WORKING_DIR_PATH)/amos/zbesh -lAMOS -l$(G2C) -lblitz -lm -o mpi_mlfma
	$(MPICC) $(INCLUDE_PATH) mpi_mlfma.o -L$(WORKING_DIR_PATH) -lMLFMA -L$(WORKING_DIR_PATH)/amos/zbesh -lAMOS -L$(WORKING_DIR_PATH)/lapack -lZGELS -l$(G2C) $(LIB_SEARCH_PATH) -lblitz -lm -o mpi_mlfma

distribute_Z_cubes: distribute_Z_cubes.o $(OBJECTS_LIBMLFMA)
	$(MPICC) $(INCLUDE_PATH) distribute_Z_cubes.o -L$(WORKING_DIR_PATH) -lMLFMA -L$(WORKING_DIR_PATH)/amos/zbesh -lAMOS -L$(WORKING_DIR_PATH)/lapack -lZGELS -l$(G2C) $(LIB_SEARCH_PATH) -lblitz -lm -o distribute_Z_cubes

mesh_functions_seb: mesh_functions_seb.o mesh.o
	$(MPICC) $(INCLUDE_PATH) mesh_functions_seb.o mesh.o readWriteBlitzArrayFromFile.o $(LIB_SEARCH_PATH) -lblitz -lm -o mesh_functions_seb

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
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) testMPI.cpp
communicateZnearBlocks.o: communicateZnearBlocks.cpp readWriteBlitzArrayFromFile.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) communicateZnearBlocks.cpp
RWGs_renumbering.o: RWGs_renumbering.cpp readWriteBlitzArrayFromFile.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) RWGs_renumbering.cpp
scatter_mesh_per_cube.o: scatter_mesh_per_cube.cpp readWriteBlitzArrayFromFile.h mesh.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) scatter_mesh_per_cube.cpp
compute_Z_near.o: compute_Z_near.cpp readWriteBlitzArrayFromFile.h Z_EJ_Z_HJ.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) compute_Z_near.cpp
mpi_mlfma.o: mpi_mlfma.cpp  Z_sparse_MLFMA.h mesh.h octtree.h octtree.cpp level.h level.cpp ./iterative/iterative.h V_E_V_H.h alpha_computation.h V_E_V_H_dipole.cpp V_E_V_H_plane.cpp EMConstants.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) mpi_mlfma.cpp
distribute_Z_cubes.o: distribute_Z_cubes.cpp octtree.h octtree.cpp level.h level.cpp EMConstants.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) distribute_Z_cubes.cpp
mesh_functions_seb.o: mesh_functions_seb.cpp mesh.cpp mesh.h readWriteBlitzArrayFromFile.cpp readWriteBlitzArrayFromFile.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) mesh_functions_seb.cpp
octtree.o: octtree.cpp interpolation.h level.h level.cpp FMM.h mesh.h octtree.h alpha_computation.h readWriteBlitzArrayFromFile.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) octtree.cpp
level.o: level.cpp interpolation.h level.h FMM.h mesh.h cube.h cube.cpp readWriteBlitzArrayFromFile.h alpha_computation.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) level.cpp
cube.o: cube.cpp FMM.h mesh.h cube.h readWriteBlitzArrayFromFile.h
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) cube.cpp
FMM.o: FMM.cpp EMConstants.h FMM.h integr_1D_X_W.h GK_triangle.h $(LIBAMOS_ZBESH_PATH)zbesh_interface.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) FMM.cpp
integr_1D_X_W.o: integr_1D_X_W.cpp GL.h 
	$(CC) $(INCLUDE_PATH) $(CFLAGS) integr_1D_X_W.cpp
interpolation.o: interpolation.cpp interpolation.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) interpolation.cpp
mesh.o: mesh.cpp mesh.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) mesh.cpp
readWriteBlitzArrayFromFile.o: readWriteBlitzArrayFromFile.cpp readWriteBlitzArrayFromFile.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) readWriteBlitzArrayFromFile.cpp
Z_EJ_Z_HJ_FS_triangles_arrays.o: Z_EJ_Z_HJ_FS_triangles_arrays.cpp EMConstants.h triangle_int.h mesh.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) Z_EJ_Z_HJ_FS_triangles_arrays.cpp
V_E_V_H_dipole.o: V_E_V_H_dipole.cpp GK_triangle.h triangle_int.h EMConstants.h V_E_V_H.h mesh.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) V_E_V_H_dipole.cpp
V_E_V_H_plane.o: V_E_V_H_plane.cpp GK_triangle.h triangle_int.h EMConstants.h V_E_V_H.h mesh.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) V_E_V_H_plane.cpp
triangle_int_FS.o: triangle_int_FS.cpp GK_triangle.h GL.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) triangle_int_FS.cpp
GK_triangle.o: GK_triangle.cpp GK_triangle.h
	$(CC) $(INCLUDE_PATH) $(CFLAGS) GK_triangle.cpp
GL.o: GL.cpp
	$(CC) $(INCLUDE_PATH) $(CFLAGS) GL.cpp



clean: 
	rm -f *.o *.a *.txt *~ *.pyc octtree mpi_mlfma testMPI communicateZnearBlocks mesh_functions_seb distribute_Z_cubes RWGs_renumbering scatter_mesh_per_cube compute_Z_near; 
	cd $(LIBAMOS_ZBESH_PATH); make clean; 
	cd $(LIBLAPACK_PATH); make clean; 
	cd $(LIBITERATIVE_PATH); make clean;

