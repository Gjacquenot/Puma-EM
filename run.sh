#!/bin/bash

DIR=$(dirname $(readlink -f $0))

n_processes=1
N_PROCESSES=${1:-$n_processes}

simu_dir="$DIR/simuDir"
SIMU_DIR=${2:-$simu_dir}

result_dir="$DIR/resultDir"
RESULT_DIR=${3:-$result_dir}

[ ! -e "${SIMU_DIR}" ] && mkdir -p ${SIMU_DIR}
rm -rf ${SIMU_DIR}/tmp*
rm -rf ${SIMU_DIR}/geo/*.msh* ${SIMU_DIR}/geo/*.txt  
mkdir -p ${SIMU_DIR}/result

echo " "
echo "Running Puma-EM on $N_PROCESSES processes. Writing in directory $SIMU_DIR."
echo " "

# load the mpi configuration variables
. ./config_mpi.sh
# custom command following mpi configuration
MPI_CMD="mpirun --hostfile $MPI_HOSTFILE -np $N_PROCESSES "

# first the mesh setup and generation. Only on one process
python code/setup_GMSH.py --simudir ${SIMU_DIR}
{ time -p ./GMSHcommand.sh; } 2> ${SIMU_DIR}/result/CPU_time_GMSH.txt

# setup of the MLFMA simulation
${MPI_CMD} python code/setup_MLFMA.py --simudir ${SIMU_DIR}

# distribution of data across processes
${MPI_CMD} ./code/MoM/distribute_Z_cubes --simudir ${SIMU_DIR}
${MPI_CMD} python code/distribute_ZChunks_and_Mesh.py --simudir ${SIMU_DIR}

# computation of the Z_near blocks
${MPI_CMD} python code/compute_Z_near_MLFMA.py --simudir ${SIMU_DIR}

# hereafter we exchange the Z_near blocks for SAI computation
# we do this in C++ as it is faster
{ time -p ${MPI_CMD} ./code/MoM/communicateZnearBlocks --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_CommunicateZnearBlocks.txt

# now computation of the SAI preconditioner
${MPI_CMD} python code/compute_SAI_precond_MLFMA.py --simudir ${SIMU_DIR}

# now the real deal: the MLFMA computation
{ time -p ${MPI_CMD} ./code/MoM/mpi_mlfma --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_MLFMA.txt

# and now the visualisation of the results
mpirun -np 1 python code/RCS_MLFMA.py --simudir ${SIMU_DIR}

# copying the result in the puma-em directory
[ ! -e "${RESULT_DIR}" ] && mkdir -p ${RESULT_DIR}
cp -r ${SIMU_DIR}/result/ ${RESULT_DIR}

# echo " "
# echo "=========================================================================="
# echo "                         SIMULATION COMPLETE! "
# echo "=========================================================================="
# echo " "
# 
# cleaning of OMPI_* variables from the environment
for i in $(env | grep OMPI_MCA |sed 's/=/ /' | awk '{print $1}')
do
	unset $i
done

exit 0
