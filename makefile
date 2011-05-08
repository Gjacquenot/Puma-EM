PROGR_DIR_PATH =  $(PWD)/code

all:
	make install_open-mpi;
install_open-mpi:
	make libs;
	killall python; killall mpi_mlfma; rm -rf tmp*;
	./run.sh;
	killall python; killall mpi_mlfma; rm -rf tmp*;
install_lam-mpi:
	make libs;
	killall python; killall mpi_mlfma; rm -rf tmp*;
	./run.sh;
	killall python; killall mpi_mlfma; rm -rf tmp*;
package:
	make clean;
	python $(PROGR_DIR_PATH)/makePackage.py;
libs: 
	cd $(PROGR_DIR_PATH)/MoM; make libs; make communicateMeshArrays; make communicateZnearBlocks; make mpi_mlfma; make mesh_functions_seb; make distribute_Z_cubes;
communicateMeshArrays:
	cd $(PROGR_DIR_PATH)/MoM; make communicateMeshArrays;
communicateZnearBlocks:
	cd $(PROGR_DIR_PATH)/MoM; make communicateZnearBlocks;
mpi_mlfma:
	cd $(PROGR_DIR_PATH)/MoM; make mpi_mlfma;
distribute_Z_cubes:
	cd $(PROGR_DIR_PATH)/MoM; make distribute_Z_cubes;
mesh_functions_seb:
	cd $(PROGR_DIR_PATH)/MoM; make mesh_functions_seb;
documentation:
	cd doc; make documentation;
clean: 
	rm -f *~ *.pyc *.txt *.out *.tar *.gz *.tgz MPIcommand.sh GMSHcommand.sh;
	cd $(PROGR_DIR_PATH); make clean; 
	cd geo; make clean;
	cd installScripts; make clean;
	cd doc; make clean;
	rm -rf Puma-EM;
	rm -rf tmp*;
	rm -rf result;
