import os

def createMPIsystemCommand(pathToProgram, programName, number_procs):
    # first analyze mpi configuration
    filename = 'config_mpi.sh'
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()
    for elem in contents:
        if "MPI_SYSTEM" in elem and "#" not in elem and "export" not in elem:
            MPI_SYSTEM = elem.split("\n")[0].split("=")[1]
        if "MPI_HOSTFILE" in elem and "#" not in elem and "export" not in elem:
            MPI_HOSTFILE = elem.split("\n")[0].split("=")[1]
    if "lam" in MPI_SYSTEM:
        MPI_command = "mpirun -np " + str(number_procs)
    elif "open-mpi" in MPI_SYSTEM and not (MPI_HOSTFILE==""):
        MPI_command = "mpirun --hostfile " + MPI_HOSTFILE + " -np " + str(number_procs)
    # now we can create the command file
    filename = 'MPIcommand.sh'
    f = open(filename, 'w')
    f.write("# automatically generated file. Cleans the OMPI_* variables from the environment\n")
    f.write("for i in $(env | grep OMPI_MCA |sed 's/=/ /' | awk '{print $1}')\n")
    f.write("do\n")
    f.write("\tunset $i\n")
    f.write("done\n")
    MPI_command += " " + os.path.join(pathToProgram, programName) + "\n"
    MPI_command += "\n"
    MPI_command += "exit 0" + "\n"
    f.write(MPI_command)
    f.close()
    os.system("chmod u+x " + filename)

def runMPIsystemCommand(pathToProgram, programName, number_procs):
    createMPIsystemCommand(pathToProgram, programName, number_procs)
    os.system("./MPIcommand.sh")
