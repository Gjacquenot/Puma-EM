#!/bin/bash

DIR=$(dirname $(readlink -f $0))

# default values
N_PROCESSES=1
INPUT_DIR="$DIR/run_in_out"
SIMU_DIR="$DIR/simuDir"

# parameters
while getopts n:i:s:h option
do
        case "${option}"
        in
                n) N_PROCESSES=$OPTARG;;
                i) INPUT_DIR=$OPTARG;;
                s) SIMU_DIR=$OPTARG;;
                h)
                        echo "Usage : ./run.sh [OPTION]..."
                        echo "  -n  Number of processes to run (default: 1)"
                        echo "  -i  Input-Output directory (default: ./run_in_out)"
                        echo "  -s  Working directory (default: ./simuDir)"
                        echo "  -h  Display this help message"
                        exit;;
        esac
done
RESULT_DIR="$INPUT_DIR/result"

echo "number of processes = $N_PROCESSES"
echo "input dir = $INPUT_DIR"
echo "simulation dir = $SIMU_DIR"
echo "result dir = $RESULT_DIR"

PYTHON_CMD="python"

[ ! -e "${SIMU_DIR}" ] && mkdir -p ${SIMU_DIR}
rm -rf ${SIMU_DIR}/tmp*
rm -rf ${SIMU_DIR}/geo/*.msh* ${SIMU_DIR}/geo/*.txt  
mkdir -p ${SIMU_DIR}/result

# JPA : debut de la partie qui sera redirigee vers le log
{

echo " "
echo "Running Puma-EM on $N_PROCESSES processes. Writing in directory $SIMU_DIR."
echo " "

# load the mpi configuration variables
. ./config_mpi.sh
# custom command following mpi configuration
MPI_CMD="mpirun --hostfile $MPI_HOSTFILE -np $N_PROCESSES "
# first the mesh setup and generation. Only on one process
${PYTHON_CMD} code/setup_GMSH.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}
{ time -p ./GMSHcommand.sh; } 2> ${SIMU_DIR}/result/CPU_time_GMSH.txt

# setup of the MLFMA simulation
${MPI_CMD} ${PYTHON_CMD} code/setup_MLFMA_folders.py --simudir ${SIMU_DIR}
${MPI_CMD} ${PYTHON_CMD} code/setup_MLFMA_excitation.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}
{ time -p ${PYTHON_CMD} code/read_MLFMA_mesh_part1.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_read_MLFMA_mesh_part1.txt
{ time -p ./code/MoM/mesh_functions_seb --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_mesh_functions_seb.txt
{ time -p ${PYTHON_CMD} code/read_MLFMA_mesh_part2.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_read_MLFMA_mesh_part2.txt
${MPI_CMD} ${PYTHON_CMD} code/setup_MLFMA_mesh.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}
${MPI_CMD} ${PYTHON_CMD} code/setup_MLFMA_poles.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}

# distribution of data across processes
{ time -p ${MPI_CMD} ./code/MoM/distribute_Z_cubes --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_distribute_Z_cubes.txt
${MPI_CMD} ${PYTHON_CMD} code/distribute_ZChunks_and_cubes.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}

# computation of the Z_near blocks
{ time -p ${MPI_CMD} ./code/MoM/compute_Z_near --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_compute_Z_near.txt
${MPI_CMD} ${PYTHON_CMD} code/compute_Z_near_MLFMA.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}

# hereafter we exchange the Z_near blocks for SAI computation
{ time -p ${MPI_CMD} ./code/MoM/communicateZnearBlocks --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_communicateZnearBlocks.txt

# now computation of the SAI preconditioner
# ${MPI_CMD} python code/compute_SAI_precond_MLFMA.py --simudir ${SIMU_DIR} --simuparams ${SIMU_PARAMS}
${MPI_CMD} ${PYTHON_CMD} code/prepare_SAI_precond_CPP.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR}
{ time -p ${MPI_CMD} ./code/MoM/compute_SAI_precond --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_compute_SAI_precond.txt

# assembling the near field interactions blocks
${MPI_CMD} ${PYTHON_CMD} code/assemble_Z_near.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR} 

# now renumbering of the RWGs for Znear and preconditioner multiplications
{ time -p ${MPI_CMD} ./code/MoM/RWGs_renumbering --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_RWGs_renumbering.txt

# now the real deal: the MLFMA computation
{ time -p ${MPI_CMD} ./code/MoM/mpi_mlfma --simudir ${SIMU_DIR}; } 2> ${SIMU_DIR}/result/CPU_time_MLFMA.txt

# and now the visualisation of the results
${PYTHON_CMD} code/RCS_MLFMA.py --inputdir ${INPUT_DIR} --simudir ${SIMU_DIR} 

# JPA : on renvoie la sortie dans un log (en plus de l'afficher dans le terminal)
} 2>&1 | tee ${SIMU_DIR}/result/output.log

# copying the result in the puma-em directory
[ ! -e "${RESULT_DIR}" ] && mkdir -p ${RESULT_DIR}
cp -r ${SIMU_DIR}/tmp0/iterative_data ${RESULT_DIR}
cp -r ${SIMU_DIR}/tmp0/octtree_data/big_cube_center_coord.txt ${RESULT_DIR}
cp -r ${SIMU_DIR}/result/* ${RESULT_DIR}
cp -r ${SIMU_DIR}/geo/*.geo ${RESULT_DIR}
rm -rf ${SIMU_DIR} GMSHcommand.sh

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

