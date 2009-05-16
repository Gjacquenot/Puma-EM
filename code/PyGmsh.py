import os, sys
import commands

def findParameter(path, targetName, parameter):
    """this function modifies a geo file through its parameters.
    If the parameter is not found it has 2 options:
    1) If TOLERANCE==True it returns an error message
    2) If TOLERANCE==False it exits from the program with an error message."""
    isGeoFileThere(path, targetName)
    fileName = os.path.join(path, targetName + '.geo')
    f = open(fileName, 'r')
    contents = f.readlines()
    f.close()
    index = 0
    for line in contents:
        if (parameter + ' = ') in line: # we try to find the line of the parameter
            IS_PARAMETER_THERE = True
            return IS_PARAMETER_THERE, index, contents
            break
        else:
            index += 1
    if index==len(contents):
        IS_PARAMETER_THERE = False
        return IS_PARAMETER_THERE, index, contents

def findParameterValue(path, targetName, parameter):
    """this function modifies a geo file through its parameters.
    If the parameter is not found it has 2 options:
    1) If TOLERANCE==True it returns an error message
    2) If TOLERANCE==False it exits from the program with an error message."""
    isGeoFileThere(path, targetName)
    IS_PARAMETER_THERE, indexParameter, contents = findParameter(path, targetName, parameter)
    if IS_PARAMETER_THERE:
        paramValueTmp = contents[indexParameter].split()[2].split(';')[0]
        return float(paramValueTmp)
    else:
        "findParameterValue(): ERROR!! Trying to find the value of nonexistent parameter", parameter, "in file", os.path.join(path, targetName + '.geo')
        sys.exit()

def write_geo(path, targetName, parameter, value):
    """this function modifies a geo file through its parameters.
    It returns an error message if the parameter is not found"""
    isGeoFileThere(path, targetName)
    IS_PARAMETER_THERE, index, contents = findParameter(path, targetName, parameter)
    if IS_PARAMETER_THERE:
        # we now replace the old value by the new value
        contents[index] = parameter + ' = ' + str(value) + ';\n'
        # we rewrite the file with the new value
        fileName = os.path.join(path, targetName + '.geo')
        f = open(fileName, 'w')
        for line in contents:
            f.write(line)
        f.close()
    else:
        print "WARNING!! Parameter", parameter, "is not in file", os.path.join(path, targetName + '.geo')

def isGeoFileThere(path, targetName):
    listOfFiles = os.listdir(path)
    geoFile = targetName + '.geo'
    for machin in listOfFiles:
        if geoFile in machin:
            return 
    print "PyGmsh: your target name", geoFile, "does not exist in", path
    print "Exiting."
    sys.exit()

def executeGmsh(path, targetName, ViewMesh):
    isGeoFileThere(path, targetName)
    #CommandString = 'gmsh -2 ' + os.path.join(path, targetName) + '.geo'
    CommandString = 'gmsh -2 -algo del2d ' + os.path.join(path, targetName) + '.geo' + ' -string "General.ExpertMode=1;"'
    print "  Meshing. Command: ", CommandString
    commands.getoutput(CommandString)
    if ViewMesh == 1:
        ViewString = 'gmsh '+ os.path.join(path, targetName) + '.msh'
        commands.getoutput(ViewString)

if __name__=="__main__":
    path = './geo'
    targetName = 'cylinder'
    TOLERANCE = True
    paramValue = findParameterValue(path, targetName, 'delta_gap', TOLERANCE)
