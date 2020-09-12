PROGR_DIR_PATH = $(PWD)/code

all:
	make install_open-mpi;
install_open-mpi:
	make libs;
	killall python; killall mpi_mlfma; rm -rf tmp*;
	./run.sh;
	killall python; killall mpi_mlfma; rm -rf tmp*;
package:
	make clean;
	python $(PROGR_DIR_PATH)/makePackage.py;
libs:
	cd $(PROGR_DIR_PATH)/MoM; make libs; \
		make communicateZnearBlocks; \
		make mpi_mlfma; \
		make mesh_functions_seb; \
		make mesh_cubes; \
		make distribute_Z_cubes; \
		make RWGs_renumbering; \
		make compute_Z_near; \
		make compute_SAI_precond;
communicateZnearBlocks:
	make - C $(PROGR_DIR_PATH)/MoM communicateZnearBlocks;
mpi_mlfma:
	make - C $(PROGR_DIR_PATH)/MoM mpi_mlfma;
distribute_Z_cubes:
	make - C $(PROGR_DIR_PATH)/MoM distribute_Z_cubes;
RWGs_renumbering:
	make - C $(PROGR_DIR_PATH)/MoM RWGs_renumbering;
compute_Z_near:
	make - C $(PROGR_DIR_PATH)/MoM compute_Z_near;
compute_SAI_precond:
	make - C $(PROGR_DIR_PATH)/MoM compute_SAI_precond;
mesh_functions_seb:
	make - C $(PROGR_DIR_PATH)/MoM mesh_functions_seb;
mesh_cubes:
	make - C $(PROGR_DIR_PATH)/MoM mesh_cubes;
documentation:
	make -C doc documentation;
clean:
	rm -rf *~ *.pyc *.txt *.out *.tar *.gz *.tgz MPIcommand.sh GMSHcommand.sh __pycache__;
	cd run_in_out; rm -rf *~ *.pyc __pycache__;
	make - C $(PROGR_DIR_PATH) clean;
	make -C geo clean;
	make -C installScripts clean;
	make -C doc clean;
	rm -rf Puma-EM;
	rm -rf tmp*;
	rm -rf result* simuDir*;
docker_run:
	if [ "${shell docker images -q pumaem 2> /dev/null}" = "" ]; then \
		docker build . -t pumaem; \
	fi
	docker run --rm -u $(shell id -u ${USER} ):$(shell id -g ${USER} ) \
       -v $(shell pwd):/opt/share -w /opt/share pumaem /bin/bash -c \
       "make CFLAGS=\"-c -O3 -fPIC -pthread -march=native -mfpmath=both\""
documentation_from_docker:
	make -C doc documentation_from_docker;
