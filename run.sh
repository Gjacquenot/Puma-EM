#!/bin/bash

echo " "
echo "Running Puma-EM..."
echo "=================="
echo " "
echo "On how many processes do you want Puma-EM to run? [default: 1]"
read N_PROCESSES

: ${N_PROCESSES:="1"}

rm -rf result

echo " "
echo "OK, Running Puma-EM on $N_PROCESSES processes..."
echo " "

# load the mpi configuration variables
. ./config_mpi.sh
# custom command following mpi configuration
if [ "$MPI_SYSTEM" = "lam-mpi" ]; then
    lamhalt
    lamboot $MPI_HOSTFILE
    COMMAND1="mpirun -np $N_PROCESSES "
else
    COMMAND1="mpirun --hostfile $MPI_HOSTFILE -np $N_PROCESSES "
fi

# cleanup. To be changed if COMPUTE_Z_NEAR is ever to have an influence again.
echo "run.sh: Erasing the tmp* directories..."
rm -rf ./tmp*
echo "Done"

# first the mesh setup and generation. Only on one process
python code/setup_GMSH.py
./GMSHcommand.sh

# setup of the MLFMA simulation
${COMMAND1} python code/setup_MLFMA.py

# distribution of data across processes
${COMMAND1} ./code/MoM/distribute_Z_cubes
${COMMAND1} python code/distribute_ZChunks_and_Mesh.py

# computation of the Z_near blocks
${COMMAND1} python code/compute_Z_near_MLFMA.py

# hereafter we exchange the Z_near blocks for SAI computation
# we do this in C++ as it is faster
${COMMAND1} ./code/MoM/communicateZnearBlocks

# now computation of the SAI preconditioner
${COMMAND1} python code/compute_SAI_precond_MLFMA.py

# now the real deal: the MLFMA computation
{ time -p ${COMMAND1} ./code/MoM/mpi_mlfma; } 2> result/CPU_time_MLFMA.txt

# and now the visualisation of the results
mpirun -np 1 python code/RCS_MLFMA.py

echo " "
echo "=========================================================================="
echo "                         SIMULATION COMPLETE! "
echo "=========================================================================="
echo " "

# cleaning of OMPI_* variables from the environment
for i in $(env | grep OMPI_MCA |sed 's/=/ /' | awk '{print $1}')
do
	unset $i
done

exit 0

