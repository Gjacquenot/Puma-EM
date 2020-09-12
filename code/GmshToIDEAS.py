import os.path, sys, math
from mesh_functions_seb import *
from PyGmsh import executeGmsh, write_geo
from EM_constants import *

def stringRightJustifiedScalarInFixedSpace(x, space):
    if isinstance(x, int):
        strTmp = str(x)
    elif isinstance(x, float):
        strTmp1 = "%.16e"%x
        strTmp = strTmp1.replace('e','D')
    blankSpace = ' ' * (space - len(strTmp))
    string = blankSpace + strTmp
    return string

def createIDEASFile(vertexes_coord, triangle_vertexes, path, targetName):
    # prologue
    f = open(os.path.join(path, targetName + '.unv'), 'w')
    f.write(stringRightJustifiedScalarInFixedSpace(-1, 6) + '\n')
    f.write(stringRightJustifiedScalarInFixedSpace(2411, 6) + '\n')
    # we first write the nodes
    V = vertexes_coord.shape[0] # the number of vertexes
    for i in range(V):
        # the node tag line
        index = int(i)
        stringToWrite = stringRightJustifiedScalarInFixedSpace(index+1, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(1, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(1, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(11, 10) + '\n'
        f.write(stringToWrite)
        # the coordinates of the node line
        stringToWrite = stringRightJustifiedScalarInFixedSpace(vertexes_coord[index, 0], 25)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(vertexes_coord[index, 1], 25)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(vertexes_coord[index, 2], 25) + '\n'
        f.write(stringToWrite)
    f.write(stringRightJustifiedScalarInFixedSpace(-1, 6) + '\n' )
    f.write(stringRightJustifiedScalarInFixedSpace(-1, 6) + '\n' )
    f.write(stringRightJustifiedScalarInFixedSpace(2412, 6) + '\n' )

    # we now write the triangles
    T = triangle_vertexes.shape[0]
    for i in range(T):
        # the triangle tag line
        index = int(i)
        stringToWrite = stringRightJustifiedScalarInFixedSpace(index+1, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(91, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(1, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(1, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(7, 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(3, 10) + '\n'
        f.write(stringToWrite)
        # the nodes indexes line
        nodesIndexes = (triangle_vertexes[index] + 1)
        stringToWrite = stringRightJustifiedScalarInFixedSpace(int(nodesIndexes[0]), 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(int(nodesIndexes[1]), 10)
        stringToWrite += stringRightJustifiedScalarInFixedSpace(int(nodesIndexes[2]), 10) + '\n'
        f.write(stringToWrite)
    f.write(stringRightJustifiedScalarInFixedSpace(-1, 6) + '\n' )


if __name__=="__main__":
    path = './geo'
    targetName = 'cubi'
    f = 2.12e9
    write_geo(path, targetName, 'lc', c/f/10)
    write_geo(path, targetName, 'lx', .1)
    write_geo(path, targetName, 'ly', .1)
    write_geo(path, targetName, 'lz', .1)
    executeGmsh(path, targetName, 0)
    z_offset = 0.0
    vertexes_coord, triangle_vertexes = read_mesh(os.path.join(path, targetName + '.msh'), z_offset)
    triangle_edgesNumbers, edgeNumber_vertexes, edgeNumber_triangles, triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = edges_computation(triangle_vertexes, vertexes_coord)
    T = len(triangle_adjacentTriangles)
    print("number of triangles T =", T)

    print("    reordering triangles for normals coherency...",)
    t0 = time.clock()
    triangles_surfaces = reorder_triangle_vertexes(triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction, triangle_vertexes, vertexes_coord)
    S = max(triangles_surfaces)+1
    time_reordering_normals = time.clock()-t0
    print("time =", time_reordering_normals, "seconds")

    t0 = time.clock()
    print("writing the mesh file in I-DEAS format...",)
    createIDEASFile(vertexes_coord, triangle_vertexes, path, targetName)
    print("time =", time.clock() - t0, "seconds")

