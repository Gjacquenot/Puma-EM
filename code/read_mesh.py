import os, sys, time
from scipy import zeros, size, compress, sort, take, put, array
from PyGmsh import executeGmsh, write_geo
#from scipy import weave
#from scipy.weave import converters

def preRead_mesh_GMSH_1(name):
    file = open(name, 'r')
    content = {}
    for line in file:
        if line[0]=='$' and 'End' not in line and 'END' not in line:
            if line.split('\n')[0] not in content:
                newFieldTmp = line.split('\n')[0]
                newField = newFieldTmp.split('$')[1]
                content[newField] = []
        elif newField in line and 'End' in line or 'END' in line:
            pass
        else:
            content[newField].append(line.split('\n')[0])
    file.close()
    return content

def read_mesh_GMSH_1(name, targetDimensions_scaling_factor, z_offset):
    """function that reads the mesh and puts it into nice arrays
       preRead does not create separate files with nodes, triangles,...
       so for big targets read_mesh_GMSH_1 can take lots of memory.

       This function is to be erased in a near future.
    """
    content = preRead_mesh_GMSH_1(name)
    if 'Nodes' in content:
        vertexes_key = 'Nodes'
    elif 'NOD' in content:
        vertexes_key = 'NOD'
    else:
        print("read_mesh.py: Error in the GMSH file format: not supported!!")
        sys.exit(1)
    V = int(content[vertexes_key][0]) # we get the number of vertexes
    vertexes_numbers = zeros(V, 'i')
    vertexes_coord = zeros( (V, 3), 'd' )

    index = 0
    for tmp in content[vertexes_key][1:]:
        tmp2 = tmp.split()
        vertexes_numbers[index] = int(tmp2[0])
        vertexes_coord[index, :] = list(map(float, tmp2[1:]))
        index += 1
    if not (targetDimensions_scaling_factor==1.0):
        vertexes_coord *= targetDimensions_scaling_factor
    vertexes_coord[:, -1] += z_offset
    
    # we check whether we have a source
    SOURCE = False
    PhysicalSurfaceNumberToName = {}
    if 'PhysicalNames' in content:
        keyPhysicalNames = 'PhysicalNames'
        NumberPhysicalNames = int(content[keyPhysicalNames][0])
        for elem in content[keyPhysicalNames][1:]:
            elem2 = elem.split()
            newKey = int(elem2[0])
            PhysicalSurfaceNumberToName[newKey] = elem2[1]
            if "source" in PhysicalSurfaceNumberToName[newKey]:
                SOURCE = True
    
    # we now extract the triangles and write them to a file
    if 'Elements' in content:
        elements_key = 'Elements'
    elif 'ELM' in content:
        elements_key = 'ELM'
    else:
        print("read_mesh.py: Error in the GMSH file format: not supported!!")
        sys.exit(1)
    N_elems = int(content[elements_key][0]) # N is the number of "elements" (border edges, facets,...)
    g = open(name + ".triangles", "w")
    #g2 = open(name + ".sourceTriangles", "w")
    T, T_src = 0, 0
    for elementTmp in content[elements_key][1:]:
        elementTmp2 = elementTmp.split()
        element = list(map(int, elementTmp2))
        if element[1]==2: # type 'planar triangle' has number 2
            listTmp = [element[3]]
            listTmp += element[-3:]
            if SOURCE and (listTmp[0] in PhysicalSurfaceNumberToName):
                T_src += 1
                for elem in listTmp:
                    #g2.write(str(elem) + " ")
                #g2.write("\n")
                    g.write(str(elem) + " ")
                g.write("\n")
            else:
                T += 1
                for elem in listTmp:
                    g.write(str(elem) + " ")
                g.write("\n")
    g.close()
    #g2.close()
    del content

    # we now read the target triangles from the name.triangles file
    g = open(name + ".triangles", 'r')
    triangles_nodes = zeros( (T + T_src, 3), 'i')
    triangles_physicalSurface = zeros(T + T_src, 'i')
    for k in range(T + T_src):
        tmp = list(map(int, g.readline().split()))
        triangles_nodes[k, :] = tmp[1:]
        triangles_physicalSurface[k] = tmp[0]
    g.close()

    vertex_number_max = max(vertexes_numbers)
    nodes_vertexes = zeros( vertex_number_max + 1, 'i' )
    put(nodes_vertexes, vertexes_numbers, range(V))
    triangles_vertexes = take(nodes_vertexes, triangles_nodes, axis=0).astype('i')
    del triangles_nodes # gain some memory, we now work only with triangles_vertexes...
    # we will now eliminate the points that we don't need
    # old non-performant memory-wise code
    #encountered_vertexesTmp = sort(array(triangles_vertexes.flat) , kind="mergesort")
    #encountered_vertexes = [encountered_vertexesTmp[0]]
    #for i in range(1, len(encountered_vertexesTmp)):
        #if (encountered_vertexesTmp[i] != encountered_vertexes[-1]):
            #encountered_vertexes.append(encountered_vertexesTmp[i])
    #max_encountered_vertexes = max(encountered_vertexes)
    # new performant memory-wise code
    max_encountered_vertexes = max(triangles_vertexes.flat)
    encountered_vertexesTmp = (zeros(max_encountered_vertexes + 1, 'i') - 1).astype('i')
    for vertex in triangles_vertexes.flat:
        encountered_vertexesTmp[vertex] = vertex
    encountered_vertexes = compress(encountered_vertexesTmp>-1, encountered_vertexesTmp, axis=0)
    del encountered_vertexesTmp
    oldVertex_to_newVertex = zeros(max_encountered_vertexes+1, 'i')
    oldVertex_to_newVertex[encountered_vertexes] = range(len(encountered_vertexes))
    triangles_vertexes = take(oldVertex_to_newVertex, triangles_vertexes, axis=0)
    vertexes_coord = take(vertexes_coord, encountered_vertexes, axis=0)
    return vertexes_coord.astype('d'), triangles_vertexes.astype('i'), triangles_physicalSurface.astype('i')

def preRead_mesh_GMSH_2(meshFile_name):
    """
       this function splits the meshFile_name.msh file into smaller entities.
       How do we do it?

       Each time we encounter a new 'entity' in the mesh file, such as 'Nodes' 
       for example, we create a new file for it, name it with the entity name,
       say 'meshFile_name.msh.Nodes', and place there all corresponding entities,
       i.e. Nodes in this example.
    """
    file = open(meshFile_name, 'r')
    content = {}
    # the special fields that can be encountered
    specialFields = ['PhysicalNames', 'Nodes', 'NOD', 'Coordinates', 'Elements', 'ELM']
    for line in file:
        # here we encounter a new entity, or field
        if line[0]=='$' and not ('End' in line or 'END' in line):
            if line.split('\n')[0] not in content:
                newFieldTmp = line.split('\n')[0]
                # we actualise newField
                newField = newFieldTmp.split('$')[1]
                content[newField] = []
                # we create the file corresponding to the entity
                fileForField = open(meshFile_name + '.' + newField, 'w')
                numberOfLines = 0
        # here we encounter the end of the new entity, or field
        elif newField in line and ('End' in line or 'END' in line):
            fileForField.close()
            content[newField] = [numberOfLines]
            pass
        else:
            if (len(content[newField])==0) and (newField in specialFields):
                content[newField] = [0]
                pass
            elif (newField in ['Elements', 'ELM']):
                elementTmp2 = line.split()
                if elementTmp2[1]=='2': # type 'planar triangle' has number 2
                  fileForField.write(line)
                  numberOfLines += 1
            else:
                fileForField.write(line)
                numberOfLines += 1
    file.close()
    return content

def read_mesh_GMSH_2(name, targetDimensions_scaling_factor, z_offset):
    """function that reads the mesh and puts it into nice arrays"""
    content = preRead_mesh_GMSH_2(name)
    print(content)
    if 'Nodes' in content:
        vertexes_key = 'Nodes'
    elif 'NOD' in content:
        vertexes_key = 'NOD'
    else:
        print("read_mesh.py: Error in the GMSH file format: not supported!!")
        sys.exit(1)
    V = int(content[vertexes_key][0]) # we get the number of vertexes
    vertexes_numbers = zeros(V, 'i')
    vertexes_coord = zeros( (V, 3), 'd' )

    index = 0
    file = open(name + '.' + vertexes_key, 'r')
    for line in file:
        tmp = line.split()
        vertexes_numbers[index] = int(tmp[0])
        vertexes_coord[index, :] = list(map(float, tmp[1:]))
        index += 1
    if not (targetDimensions_scaling_factor==1.0):
        vertexes_coord *= targetDimensions_scaling_factor
    vertexes_coord[:, -1] += z_offset
    file.close()

    # we check whether we have a source
    SOURCE = False
    PhysicalSurfaceNumberToName = {}
    if 'PhysicalNames' in content:
        keyPhysicalNames = 'PhysicalNames'
        NumberPhysicalNames = content[keyPhysicalNames][0]
        file = open(name + '.' + keyPhysicalNames, 'r')
        for line in file:
            elem2 = line.split()
            newKey = int(elem2[0])
            PhysicalSurfaceNumberToName[newKey] = elem2[1]
            if "source" in PhysicalSurfaceNumberToName[newKey]:
                SOURCE = True
        file.close()

    # we now extract the triangles and write them to a file
    if 'Elements' in content:
        elements_key = 'Elements'
    elif 'ELM' in content:
        elements_key = 'ELM'
    else:
        print("read_mesh.py: Error in the GMSH file format: not supported!!")
        sys.exit(1)
    T = content[elements_key][0] # N is the number of triangles
    del content

    g = open(name + "." + elements_key, 'r')
    triangles_nodes = zeros( (T, 3), 'i')
    triangles_physicalSurface = zeros(T, 'i')
    index = 0
    for line in g:
        tmp = list(map(int, line.split()))
        triangles_nodes[index, :] = tmp[-3:]
        triangles_physicalSurface[index] = tmp[3]
        index += 1
    g.close()
    
    # indexes of elements in Python/C++ arrays start at 0.
    # However, triangles_nodes don't necessarily start at 0.
    # So the following 4 lines correct that.
    vertex_number_max = max(vertexes_numbers)
    nodes_vertexes = zeros( vertex_number_max + 1, 'i' )
    put(nodes_vertexes, vertexes_numbers, range(V))
    triangles_vertexes = take(nodes_vertexes, triangles_nodes, axis=0).astype('i')
    # we will now eliminate the points that we don't need
#    del triangles_nodes # gain some memory, we now work only with triangles_vertexes...
#    max_encountered_vertexes = max(triangles_vertexes.flat)
#    encountered_vertexesTmp = (zeros(max_encountered_vertexes + 1, 'i') - 1).astype('i')
#    wrapping_code = """
#    for (int i=0 ; i<triangles_vertexes.extent(0) ; i++) {
#      for (int j=0 ; j<triangles_vertexes.extent(1) ; j++) {
#        const int vertex = triangles_vertexes(i, j);
#        encountered_vertexesTmp(vertex) = vertex;
#      }
#    }
#    """
#    weave.inline(wrapping_code,
#                 ['encountered_vertexesTmp', 'triangles_vertexes'],
#                 type_converters = converters.blitz,
#                 include_dirs = [],
#                 library_dirs = [],
#                 libraries = [],
#                 headers = ['<iostream>'],
#                 compiler = 'gcc',
#                 extra_compile_args = ['-O3', '-pthread', '-w'])
#    encountered_vertexes = compress(encountered_vertexesTmp>-1, encountered_vertexesTmp, axis=0)
#    del encountered_vertexesTmp
#    oldVertex_to_newVertex = zeros(max_encountered_vertexes+1, 'i')
#    oldVertex_to_newVertex[encountered_vertexes] = range(len(encountered_vertexes))
#    triangles_vertexes = take(oldVertex_to_newVertex, triangles_vertexes, axis=0)
#    vertexes_coord = take(vertexes_coord, encountered_vertexes, axis=0)
    return vertexes_coord.astype('d'), triangles_vertexes.astype('i'), triangles_physicalSurface.astype('i')

def preRead_mesh_GiD(meshFile_name):
    """
       GiD is a mesh generator that you can google on the Internet.

       This function splits the meshFile_name.msh file into smaller entities.
       How do we do it?

       Each time we encounter a new 'entity' in the mesh file, such as 'Nodes' 
       for example, we create a new file for it, name it with the entity name,
       say 'meshFile_name.msh.Nodes', and place there all corresponding entities,
       i.e. Nodes in this example.
    """
    file = open(meshFile_name, 'r')
    file.readline() # skip the first line
    file.readline() # skip the Coordinates line
    content = {}
    coordinatesFile = open(meshFile_name + '.' + 'Nodes', 'w')
    numberOfLines = 0
    for line in file:
        if 'end' in line:
            break
        else:
            numberOfLines += 1
            coordinatesFile.write(line)
    content['Nodes'] = [numberOfLines]
    coordinatesFile.close()
    # we now read and rewrite the 'Elements'
    elementsFile = open(meshFile_name + '.' + 'Elements', 'w')
    numberOfLines = 0
    for lineTmp in file:
        line = lineTmp.split()
        if 'Elements' in line or len(line)==0:
            pass
        elif 'end' in lineTmp:
            break
        else:
            numberOfLines += 1
            elementsFile.write(lineTmp)
    content['Elements'] = [numberOfLines]
    elementsFile.close()
    # final moves
    file.close()
    return content

def read_mesh_GiD(name, targetDimensions_scaling_factor, z_offset):
    """function that reads the mesh and puts it into nice arrays"""
    content = preRead_mesh_GiD(name)
    print(content)
    vertexes_key = 'Nodes'

    V = int(content[vertexes_key][0]) # we get the number of vertexes
    vertexes_numbers = zeros(V, 'i')
    vertexes_coord = zeros( (V, 3), 'd' )

    index = 0
    file = open(name + '.' + vertexes_key, 'r')
    for line in file:
        tmp = line.split()
        vertexes_numbers[index] = int(tmp[0])
        vertexes_coord[index, :] = list(map(float, tmp[1:]))
        index += 1
    if not (targetDimensions_scaling_factor==1.0):
        vertexes_coord *= targetDimensions_scaling_factor
    vertexes_coord[:, -1] += z_offset
    file.close()

    # we now extract the triangles and write them to a file
    if 'Elements' in content:
        elements_key = 'Elements'
    else:
        print("read_mesh.py: Error in the GiD file format: not supported!!")
        sys.exit(1)
    T = content[elements_key][0] # N is the number of triangles
    del content

    g = open(name + "." + elements_key, 'r')
    triangles_nodes = zeros( (T, 3), 'i')
    triangles_physicalSurface = zeros(T, 'i')
    index = 0
    for line in g:
        tmp = list(map(int, line.split()))
        triangles_nodes[index, :] = tmp[-3:]
        #triangles_physicalSurface[index] = tmp[3]
        triangles_physicalSurface[index] = 0 # not like GMSH here: GMSH triangles have a physical surface number
        index += 1
    g.close()
    
    # indexes of elements in Python/C++ arrays start at 0.
    # However, triangles_nodes don't necessarily start at 0.
    # So the following 4 lines correct that.
    vertex_number_max = max(vertexes_numbers)
    nodes_vertexes = zeros( vertex_number_max + 1, 'i' )
    put(nodes_vertexes, vertexes_numbers, range(V))
    triangles_vertexes = take(nodes_vertexes, triangles_nodes, axis=0).astype('i')
    # we will now eliminate the points that we don't need
#    del triangles_nodes # gain some memory, we now work only with triangles_vertexes...
#    max_encountered_vertexes = max(triangles_vertexes.flat)
#    encountered_vertexesTmp = (zeros(max_encountered_vertexes + 1, 'i') - 1).astype('i')
#    wrapping_code = """
#    for (int i=0 ; i<triangles_vertexes.extent(0) ; i++) {
#      for (int j=0 ; j<triangles_vertexes.extent(1) ; j++) {
#        const int vertex = triangles_vertexes(i, j);
#        encountered_vertexesTmp(vertex) = vertex;
#      }
#    }
#    """
#    weave.inline(wrapping_code,
#                 ['encountered_vertexesTmp', 'triangles_vertexes'],
#                 type_converters = converters.blitz,
#                 include_dirs = [],
#                 library_dirs = [],
#                 libraries = [],
#                 headers = ['<iostream>'],
#                 compiler = 'gcc',
#                 extra_compile_args = ['-O3', '-pthread', '-w'])
#    encountered_vertexes = compress(encountered_vertexesTmp>-1, encountered_vertexesTmp, axis=0)
#    del encountered_vertexesTmp
#    oldVertex_to_newVertex = zeros(max_encountered_vertexes+1, 'i')
#    oldVertex_to_newVertex[encountered_vertexes] = range(len(encountered_vertexes))
#    triangles_vertexes = take(oldVertex_to_newVertex, triangles_vertexes, axis=0)
#    vertexes_coord = take(vertexes_coord, encountered_vertexes, axis=0)
    return vertexes_coord.astype('d'), triangles_vertexes.astype('i'), triangles_physicalSurface.astype('i')

def preRead_mesh_ANSYS(meshFile_path):
    """
       ANSYS is a multiphysics program that can mesh geometries.

       This function splits the meshFile_name.msh file into smaller entities.
       How do we do it?

       Each time we encounter a new 'entity' in the mesh file, such as 'Nodes' 
       for example, we create a new file for it, name it with the entity name,
       say 'meshFile_name.msh.Nodes', and place there all corresponding entities,
       i.e. Nodes in this example.
    """
    # we first open the nodes coordinates file from which we will read
    # NLIST.lis is the default ANSYS name.
    fileToRead = open(os.path.join(meshFile_path, 'NLIST.lis'), 'r')
    # we also open the file into which we will write the reformatted nodes
    coordinatesFile = open(os.path.join(meshFile_path, 'NLIST.lis.Nodes'), 'w')
    # various initialisations
    content = {}
    numberOfLines = 0
    # and now the read and write operations
    for currentLine in fileToRead:
        words = currentLine.split()
        # we pass the line if it is empty. If non empty, we test for numeral entries.
        if len(words) > 0:
            word = words[0]
            if (word[0]>='0' and word[0]<='9'): # we then have an integer
                numberOfLines += 1
                lineToWrite = words[0] + ' ' + words[1] + ' ' + words[2] + ' ' + words[3] + '\n'
                coordinatesFile.write(lineToWrite)
    content['Nodes'] = [numberOfLines]
    fileToRead.close()
    coordinatesFile.close()
    # we now read the 'Elements'
    # ELIST.lis is the default ANSYS name.
    fileToRead = open(os.path.join(meshFile_path, 'ELIST.lis'), 'r')
    elementsFile = open(os.path.join(meshFile_path, 'ELIST.lis.Elements'), 'w')
    numberOfLines = 0
    for currentLine in fileToRead:
        words = currentLine.split()
        # we pass the line if it is empty. If non empty, we test for numeral entries.
        if len(words) > 0:
            word = words[0]
            if (word[0]>='0' and word[0]<='9'): # we then have an integer
                numberOfLines += 1
                lineToWrite = words[0] + ' ' + words[6] + ' ' + words[7] + ' ' + words[8] + '\n'
                elementsFile.write(lineToWrite)
    content['Elements'] = [numberOfLines]
    elementsFile.close()
    fileToRead.close()
    return content


def read_mesh_ANSYS(path, name, targetDimensions_scaling_factor, z_offset):
    """function that reads the mesh and puts it into nice arrays"""
    content = preRead_mesh_ANSYS(path) # name doesnt count here, since ANSYS has default names
    print(content)
    vertexes_key = 'Nodes'

    V = int(content[vertexes_key][0]) # we get the number of vertexes
    vertexes_numbers = zeros(V, 'i')
    vertexes_coord = zeros( (V, 3), 'd' )

    index = 0
    fileToRead = open(os.path.join(path, 'NLIST.lis' + '.' + vertexes_key), 'r')
    for line in fileToRead:
        tmp = line.split()
        vertexes_numbers[index] = int(tmp[0])
        vertexes_coord[index, :] = list(map(float, tmp[1:]))
        index += 1
    if not (targetDimensions_scaling_factor==1.0):
        vertexes_coord *= targetDimensions_scaling_factor
    vertexes_coord[:, -1] += z_offset
    fileToRead.close()

    # we now extract the triangles and write them to a file
    if 'Elements' in content:
        elements_key = 'Elements'
    else:
        print("read_mesh.py: Error in the ANSYS file format: not supported!!")
        sys.exit(1)
    T = content[elements_key][0] # N is the number of triangles
    del content

    g = open(os.path.join(path, 'ELIST.lis' + "." + elements_key), 'r')
    triangles_nodes = zeros( (T, 3), 'i')
    triangles_physicalSurface = zeros(T, 'i')
    index = 0
    for line in g:
        tmp = list(map(int, line.split()))
        triangles_nodes[index, :] = tmp[-3:]
        #triangles_physicalSurface[index] = tmp[3]
        triangles_physicalSurface[index] = 0 # not like GMSH here: GMSH triangles have a physical surface number
        index += 1
    g.close()
    
    # indexes of elements in Python/C++ arrays start at 0.
    # However, triangles_nodes don't necessarily start at 0.
    # So the following 4 lines correct that.
    vertex_number_max = max(vertexes_numbers)
    nodes_vertexes = zeros( vertex_number_max + 1, 'i' )
    put(nodes_vertexes, vertexes_numbers, range(V))
    triangles_vertexes = take(nodes_vertexes, triangles_nodes, axis=0).astype('i')
    # we will now eliminate the points that we don't need
#    del triangles_nodes # gain some memory, we now work only with triangles_vertexes...
#    max_encountered_vertexes = max(triangles_vertexes.flat)
#    encountered_vertexesTmp = (zeros(max_encountered_vertexes + 1, 'i') - 1).astype('i')
#    wrapping_code = """
#    for (int i=0 ; i<triangles_vertexes.extent(0) ; i++) {
#      for (int j=0 ; j<triangles_vertexes.extent(1) ; j++) {
#        const int vertex = triangles_vertexes(i, j);
#        encountered_vertexesTmp(vertex) = vertex;
#      }
#    }
#    """
#    weave.inline(wrapping_code,
#                 ['encountered_vertexesTmp', 'triangles_vertexes'],
#                 type_converters = converters.blitz,
#                 include_dirs = [],
#                 library_dirs = [],
#                 libraries = [],
#                 headers = ['<iostream>'],
#                 compiler = 'gcc',
#                 extra_compile_args = ['-O3', '-pthread', '-w'])
#    encountered_vertexes = compress(encountered_vertexesTmp>-1, encountered_vertexesTmp, axis=0)
#    del encountered_vertexesTmp
#    oldVertex_to_newVertex = zeros(max_encountered_vertexes+1, 'i')
#    oldVertex_to_newVertex[encountered_vertexes] = range(len(encountered_vertexes))
#    triangles_vertexes = take(oldVertex_to_newVertex, triangles_vertexes, axis=0)
#    vertexes_coord = take(vertexes_coord, encountered_vertexes, axis=0)
    return vertexes_coord.astype('d'), triangles_vertexes.astype('i'), triangles_physicalSurface.astype('i')



if __name__=="__main__":
    path = './geo'
    targetName = 'sphere'
    write_geo(path, targetName, 'lc', 0.05)
    write_geo(path, targetName, 'lx', 0.051)
    write_geo(path, targetName, 'ly', 0.1)
    write_geo(path, targetName, 'lz', 0.2)
    write_geo(path, targetName, 'frill_width', 0.02)
    executeGmsh(path, targetName, 0)
    z_offset = 0.0
    targetDimensions_scaling_factor = 1.0
    t0 = time.time()
    vertexes_coord_1, triangles_vertexes_1, triangles_physicalSurface_1 = read_mesh_GMSH_1(os.path.join(path, targetName) + '.msh', targetDimensions_scaling_factor, z_offset)
    print("time for classical *.msh file reading = " + str(time.time() - t0))
    t0 = time.time()
    vertexes_coord_2, triangles_vertexes_2, triangles_physicalSurface_2 = read_mesh_GMSH_2(os.path.join(path, targetName) + '.msh', targetDimensions_scaling_factor, z_offset)
    print("time for new *.msh file reading = " + str(time.time() - t0))
    
    print
    print("difference between python and C++ code. If results different than 0, there is a problem.")
    print(str(sum(abs(vertexes_coord_1 - vertexes_coord_2))))
    print(str(sum(abs(triangles_vertexes_1 - triangles_vertexes_2))))
    print(str(sum(abs(triangles_physicalSurface_1 - triangles_physicalSurface_2))))
    
    vertexes_coord_GiD, triangles_vertexes_GiD, triangles_physicalSurface_GiD = read_mesh_GiD(os.path.join(path, 'aaa1') + '.msh', targetDimensions_scaling_factor, z_offset)

    content = preRead_mesh_ANSYS('./geo')
    vertexes_coord, triangles_vertexes, triangles_physicalSurface = read_mesh_ANSYS(path, "whatever", targetDimensions_scaling_factor, z_offset)
    print(vertexes_coord.shape, triangles_vertexes.shape, triangles_physicalSurface.shape)

