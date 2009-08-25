import os.path, sys, time
from scipy import zeros, ones, arange, array, take, reshape, sort, argsort, put, sum, compress, nonzero, prod
from scipy import mean, sqrt
from read_mesh import read_mesh_GMSH_1, read_mesh_GMSH_2, read_mesh_GIS
from Cubes import cube_lower_coord_computation, RWGNumber_cubeNumber_computation, cubeIndex_RWGNumbers_computation, findCubeNeighbors, cubes_indexes_to_numbers_computation, write_cubes
from PyGmsh import executeGmsh, write_geo, findParameter, findParameterValue
from ReadWriteBlitzArray import readBlitzArrayFromDisk, read1DBlitzArrayFromDisk, readIntFromDisk, writeBlitzArrayToDisk, writeScalarToDisk
from EM_constants import *
from mesh_functions_seb import *
from mesh_functions_seb_C import *
import copy

class MeshClass:
    def __init__(self, path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction, meshFormat, meshFileTermination):
	
	self.path, self.name, self.geoName = path, targetName + meshFileTermination, targetName + '.geo'
        self.targetName = targetName
        self.z_offset = z_offset
        self.targetDimensions_scaling_factor = targetDimensions_scaling_factor
        self.languageForMeshConstruction = languageForMeshConstruction
        self.meshFormat = meshFormat
        self.meshFileTermination = meshFileTermination

    def constructFromGmshFile(self):
        # we check to see if there is a delta_gap parameter in the geo file
        if self.meshFormat == 'GMSH':
            self.DELTA_GAP, tmp1, tmp2 = findParameter(self.path, self.targetName, "delta_gap")
            if self.DELTA_GAP:
                ## a delta gap should always be defined in the *.geo file as being located
                ## at the origin r = [0,0,0], and directed following z"""
                self.delta_gap = findParameterValue(self.path, self.targetName, "delta_gap")
                print "There is a delta gap in file", self.geoName
                print "The size of the delta gap is", self.delta_gap, "m"
        else:
            # else we don't look for a delta gap in the file
            pass
        print "  Python construction of the Mesh object"
        sys.stdout.flush()

        print "  reading",  os.path.join(self.path, self.name), "...",
        t0 = time.clock()
        if self.meshFormat == 'GMSH':
            #self.vertexes_coord, self.triangle_vertexes, self.triangles_physicalSurface = read_mesh_GMSH_1(os.path.join(self.path, self.name), self.targetDimensions_scaling_factor, self.z_offset)
            self.vertexes_coord, self.triangle_vertexes, self.triangles_physicalSurface = read_mesh_GMSH_2(os.path.join(self.path, self.name), self.targetDimensions_scaling_factor, self.z_offset)
        elif self.meshFormat == 'GIS':
            self.vertexes_coord, self.triangle_vertexes, self.triangles_physicalSurface = read_mesh_GIS(os.path.join(self.path, self.name), self.targetDimensions_scaling_factor, self.z_offset)
        else:
            print "meshClass.py : error on the mesh format. Enter a correct one please."
        self.time_reading = time.clock()-t0
        print "reading mesh time =", self.time_reading, "seconds"
        self.T = self.triangle_vertexes.shape[0]
        print "  number of triangles =", self.T
        print "  edges classification..."
        sys.stdout.flush()
        
        if self.languageForMeshConstruction=="C" or self.languageForMeshConstruction=="C++":
            t0 = time.clock()
            self.triangles_surfaces, self.IS_CLOSED_SURFACE, self.RWGNumber_signedTriangles, self.RWGNumber_edgeVertexes = edges_computation_C(self.triangle_vertexes, self.vertexes_coord)
            self.N_RWG = self.RWGNumber_signedTriangles.shape[0]
            self.S = len(self.IS_CLOSED_SURFACE)
            print "  test of the closed surfaces :", self.IS_CLOSED_SURFACE
            self.RWGNumber_oppVertexes = RWGNumber_oppVertexes_computation_C(self.RWGNumber_signedTriangles, self.RWGNumber_edgeVertexes, self.triangle_vertexes)
            self.time_edges_classification = time.clock()-t0
            print "  edges classification cumulated time =", self.time_edges_classification, "seconds"
            self.time_effective_RWG_functions_computation = self.time_edges_classification
            print "    Number of RWG =", self.N_RWG

        else:
            t0 = time.clock()
            edgeNumber_vertexes, edgeNumber_triangles, triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = edges_computation(self.triangle_vertexes, self.vertexes_coord)
            self.time_edges_classification = time.clock()-t0
            print "  edges classification cumulated time =", self.time_edges_classification, "seconds"

            print "  reordering triangles for normals coherency..."
            sys.stdout.flush()
            t0 = time.clock()
            self.triangles_surfaces = reorder_triangle_vertexes(triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction, self.triangle_vertexes, self.vertexes_coord)
            self.S = max(self.triangles_surfaces)+1
            self.time_reordering_normals = time.clock()-t0
            print "  cumulated time =", self.time_reordering_normals, "seconds"
    
            print "  checking the closed and open surfaces...",
            sys.stdout.flush()
            t0 = time.clock()
            self.IS_CLOSED_SURFACE, self.connected_surfaces, self.potential_closed_surfaces = is_surface_closed(self.triangles_surfaces, edgeNumber_triangles)
            print "  test of the closed surfaces :", self.IS_CLOSED_SURFACE
            print "  connected surfaces : ", self.connected_surfaces
            print "  potential closed surfaces : ", self.potential_closed_surfaces
        
            print "  computing the effective RWG functions and their opposite vertexes..."
            sys.stdout.flush()
            t0 = time.clock()
            self.RWGNumber_signedTriangles, self.RWGNumber_edgeVertexes, self.N_edges, self.N_RWG = RWGNumber_signedTriangles_computation(edgeNumber_triangles, edgeNumber_vertexes, self.triangles_surfaces, self.IS_CLOSED_SURFACE, self.triangle_vertexes, self.vertexes_coord)
            del edgeNumber_vertexes
            self.RWGNumber_oppVertexes = RWGNumber_oppVertexes_computation(self.RWGNumber_signedTriangles, self.RWGNumber_edgeVertexes, self.triangle_vertexes)
            # memory-economic way for computing average_RWG_length
            self.time_effective_RWG_functions_computation =  time.clock() - t0
            print "  effective RWG functions computation cumulated time =", self.time_effective_RWG_functions_computation
            print "  Number of edges =", self.N_edges
            print "  Number of RWG =", self.N_RWG
            sys.stdout.flush()
        
        self.compute_RWG_CFIE_OK()
        if self.N_RWG<1e4:
            stride = 1
        else:
            stride = self.N_RWG/100
        self.average_RWG_length = compute_RWG_meanEdgeLength(self.vertexes_coord, self.RWGNumber_edgeVertexes, stride)
        
    def normals_check(self):
        triangles_areas, triangles_normals = triangles_areas_normals_computation(self.vertexes_coord, self.triangle_vertexes, self.triangles_surfaces)
        self.test_normal_integral = zeros(self.S, 'd')
        for s in range(self.S):
            does_t_belongs_to_s = (self.triangles_surfaces == s)
            self.test_normal_integral[s] = sum(sum(compress(does_t_belongs_to_s, triangles_areas, 0) * compress(does_t_belongs_to_s, triangles_normals, 0), 1))

    def compute_RWG_CFIE_OK(self):
        RWGNumber_CFIE_OK_tmp1 = take(self.triangles_surfaces, self.RWGNumber_signedTriangles, axis=0)
        RWGNumber_CFIE_OK_tmp2 = take(self.IS_CLOSED_SURFACE, RWGNumber_CFIE_OK_tmp1, axis=0)
        # following the Taskinen et al. paper in PIER
        # we can have a CFIE on a junction straddling a dielectric and metallic body
        self.RWGNumber_CFIE_OK = ((sum(RWGNumber_CFIE_OK_tmp2, axis=1)>=1) * 1).astype('i')
        # We cannot have M on a junction between a dielectric and metallic body!
        # The following expression is lacking the fact that a surface can be metallic
        # or dielectric. If metallic, then there is no M current, even if surface is closed
        self.RWGNumber_M_CURRENT_OK = ((sum(RWGNumber_CFIE_OK_tmp2, axis=1)>1) * 1).astype('i')
    
    def write_normals(self):
        print "writing normals to a file"
        name = 'normals.pos'
        print "triangles_areas and triangles_normals computation...",
        t0 = time.clock()
        triangles_centroids = triangles_centroids_computation(self.vertexes_coord, self.triangle_vertexes)
        triangles_areas, triangles_normals = triangles_areas_normals_computation(self.vertexes_coord, self.triangle_vertexes, self.triangles_surfaces)
        self.time_triangles_normals_comp = time.clock()-t0
        print "time =", self.time_triangles_normals_comp, "seconds"
        write_normals(os.path.join(self.path, name), triangles_centroids, triangles_normals, self.triangles_surfaces, -1)

    def cubes_data_computation(self, a):
        #print "Entering cubes_data_computation............."
        self.a = a
        self.max_N_cubes_1D, self.N_levels, self.big_cube_lower_coord, self.big_cube_center_coord = cube_lower_coord_computation(a, self.vertexes_coord)
        self.N_levels = max(self.N_levels, 2)
        #print "N levels =", self.N_levels, ", max N cubes 1D =", self.max_N_cubes_1D
        RWGNumber_edgeCentroidCoord = compute_RWGNumber_edgeCentroidCoord(self.vertexes_coord, self.RWGNumber_edgeVertexes)
        RWGNumber_cube, RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord = RWGNumber_cubeNumber_computation(a, self.max_N_cubes_1D, self.big_cube_lower_coord, RWGNumber_edgeCentroidCoord)
        self.cubes_RWGsNumbers, self.cubes_lists_RWGsNumbers, self.cubes_Ei, self.cubes_centroids = cubeIndex_RWGNumbers_computation(RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord)
        self.C = self.cubes_centroids.shape[0]
        self.cubesNeighborsIndexes = findCubeNeighbors(self.max_N_cubes_1D, self.big_cube_lower_coord, self.cubes_centroids, self.a, self.N_levels)
        self.cubes_indexes_to_numbers = cubes_indexes_to_numbers_computation(self.a, self.big_cube_lower_coord, self.cubes_centroids, self.N_levels)
        #print "Average number of RWGs per cube:", mean(self.cubes_Ei)
        #print "Exiting cubes_data_computation.............."

    def saveToDisk(self, path):
        """this function writes to disk the arrays necessary for MLFMA to work."""
        writeBlitzArrayToDisk(self.cubes_centroids, os.path.join(path, 'cubes_centroids') + '.txt')
        writeBlitzArrayToDisk(self.cubes_RWGsNumbers, os.path.join(path, 'cubes_RWGsNumbers') + '.txt')
        writeBlitzArrayToDisk(self.vertexes_coord, os.path.join(path, 'vertexes_coord') + '.txt')
        writeBlitzArrayToDisk(self.triangle_vertexes, os.path.join(path, 'triangle_vertexes') + '.txt')
        writeBlitzArrayToDisk(self.triangles_surfaces, os.path.join(path, 'triangles_surfaces') + '.txt')
        writeBlitzArrayToDisk(self.IS_CLOSED_SURFACE, os.path.join(path, 'isClosedSurface') + '.txt')
        writeBlitzArrayToDisk(self.RWGNumber_CFIE_OK, os.path.join(path, 'RWGNumber_CFIE_OK') + '.txt')
        writeBlitzArrayToDisk(self.RWGNumber_M_CURRENT_OK, os.path.join(path, 'RWGNumber_M_CURRENT_OK') + '.txt')
        writeBlitzArrayToDisk(self.RWGNumber_signedTriangles, os.path.join(path, 'RWGNumber_signedTriangles') + '.txt')
        writeBlitzArrayToDisk(self.RWGNumber_edgeVertexes, os.path.join(path, 'RWGNumber_edgeVertexes') + '.txt')
        writeBlitzArrayToDisk(self.RWGNumber_oppVertexes, os.path.join(path, 'RWGNumber_oppVertexes') + '.txt')
        writeBlitzArrayToDisk(self.big_cube_lower_coord, os.path.join(path, 'big_cube_lower_coord') + '.txt')
        writeBlitzArrayToDisk(self.big_cube_center_coord, os.path.join(path, 'big_cube_center_coord') + '.txt')
        # now we write the scalar values
        writeScalarToDisk(self.C, os.path.join(path, "C.txt"))
        writeScalarToDisk(self.N_RWG, os.path.join(path, "E.txt"))
        writeScalarToDisk(self.T, os.path.join(path, "T.txt"))
        writeScalarToDisk(self.S, os.path.join(path, "S.txt"))
        writeScalarToDisk(self.vertexes_coord.shape[0], os.path.join(path, "V.txt"))
        writeScalarToDisk(self.N_levels, os.path.join(path, "N_levels.txt"))
        writeScalarToDisk(self.cubes_RWGsNumbers.shape[1], os.path.join(path, "MAX_N_RWG_per_cube.txt"))

    def constructFromSavedArrays(self, path):
        self.C = readIntFromDisk(os.path.join(path, "C.txt"))
        self.E = readIntFromDisk(os.path.join(path, "E.txt"))
        self.S = readIntFromDisk(os.path.join(path, "S.txt"))
        self.T = readIntFromDisk(os.path.join(path, "T.txt"))
        self.V = readIntFromDisk(os.path.join(path, "V.txt"))
        self.N_levels = readIntFromDisk(os.path.join(path, "N_levels.txt"))
        self.MAX_N_RWG_per_cube = readIntFromDisk(os.path.join(path, "MAX_N_RWG_per_cube.txt"))
        self.N_RWG = self.E
        # now the arrays
        self.cubes_centroids = readBlitzArrayFromDisk(os.path.join(path, "cubes_centroids.txt"), self.C, 3, 'd')
        self.cubes_RWGsNumbers = readBlitzArrayFromDisk(os.path.join(path, "cubes_RWGsNumbers.txt"), self.C, self.MAX_N_RWG_per_cube, 'i')
        self.vertexes_coord = readBlitzArrayFromDisk(os.path.join(path, "vertexes_coord.txt"), self.V, 3, 'd')
        self.triangle_vertexes = readBlitzArrayFromDisk(os.path.join(path, "triangle_vertexes.txt"), self.T, 3, 'i')
        self.triangles_surfaces = read1DBlitzArrayFromDisk(os.path.join(path, "triangles_surfaces.txt"), 'i')
        self.RWGNumber_CFIE_OK = read1DBlitzArrayFromDisk(os.path.join(path, "RWGNumber_CFIE_OK.txt"), 'i')
        self.RWGNumber_M_CURRENT_OK = read1DBlitzArrayFromDisk(os.path.join(path, "RWGNumber_M_CURRENT_OK.txt"), 'i')
        self.RWGNumber_signedTriangles = readBlitzArrayFromDisk(os.path.join(path, "RWGNumber_signedTriangles.txt"), self.E, 2, 'i')
        self.RWGNumber_edgeVertexes = readBlitzArrayFromDisk(os.path.join(path, "RWGNumber_edgeVertexes.txt"), self.E, 2, 'i')
        self.RWGNumber_oppVertexes = readBlitzArrayFromDisk(os.path.join(path, "RWGNumber_oppVertexes.txt"), self.E, 2, 'i')
        self.big_cube_lower_coord = read1DBlitzArrayFromDisk(os.path.join(path, "big_cube_lower_coord.txt"), 'd')
        self.big_cube_center_coord = read1DBlitzArrayFromDisk(os.path.join(path, "big_cube_center_coord.txt"), 'd')
        self.average_RWG_length = mean(compute_RWGNumber_edgeLength(self.vertexes_coord, self.RWGNumber_edgeVertexes))
        # checking the normals
        self.IS_CLOSED_SURFACE = read1DBlitzArrayFromDisk(os.path.join(path, "isClosedSurface.txt"), 'i')

    def sizeOfZnear_per_Cube(self, cubeNumber):
        test_edges_numbers = self.cubes_lists_edges_numbers[int(cubeNumber)]
        Number_of_src_edges = 0
        for j in self.cubesNeighborsIndexes[int(cubeNumber)]:
            Number_of_src_edges += len(self.cubes_lists_edges_numbers[int(j)])
        return len(test_edges_numbers), Number_of_src_edges

    def list_testSrc_RWGs_computation(self, cubeNumber):
        # we now construct the source and observer edges lists
        list_of_test_edges_numbers = self.cubes_lists_RWGsNumbers[cubeNumber].tolist()
        list_of_src_edges_numbers = []
        for j in self.cubesNeighborsIndexes[cubeNumber]:
            list_of_src_edges_numbers += self.cubes_lists_RWGsNumbers[int(j)].tolist()
        return list_of_test_edges_numbers, list_of_src_edges_numbers

    def computeCubeLocalArrays(self, cubeNumber):
        """this function writes to disk the arrays pertaining to a cube.
        this should allow for localized near field and preconditioner computation, 
        i.e. only with local mesh info."""
        list_of_test_RWGsNumbers, list_of_src_RWGsNumbers = self.list_testSrc_RWGs_computation(cubeNumber)
        N_RWG_test, N_RWG_src = len(list_of_test_RWGsNumbers), len(list_of_src_RWGsNumbers)
        # computation of the local RWGs and their characteristics
        localTestSrcRWGNumber_signedTriangles = take(self.RWGNumber_signedTriangles, list_of_src_RWGsNumbers, axis=0)
        # we have to know for which basis function we can apply the CFIE. It is ok if at least one
        # triangle of the RWG is on a closed surface.
        localTestRWGNumber_CFIE_OK = take(self.RWGNumber_CFIE_OK, list_of_test_RWGsNumbers, axis=0).astype('i')
        localSrcRWGNumber_M_CURRENT_OK = take(self.RWGNumber_M_CURRENT_OK, list_of_src_RWGsNumbers, axis=0).astype('i')
        # RWGNumber_vertexesCoord
        localTestSrcRWGNumber_vertexesTmp = take(self.RWGNumber_edgeVertexes, list_of_src_RWGsNumbers, axis=0)
        #localTestSrcRWGNumber_vertexesCoord = take(self.vertexes_coord, localTestSrcRWGNumber_vertexes, axis=0)
        # RWGNumber_oppVertexesCoord
        localTestSrcRWGNumber_oppVertexesTmp = take(self.RWGNumber_oppVertexes, list_of_src_RWGsNumbers, axis=0)
        #localTestSrcRWGNumber_oppVertexesCoord = take(self.vertexes_coord, localTestSrcRWGNumber_oppVertexes, axis=0)

        localTestSrcRWGNumber_vertexes = zeros((localTestSrcRWGNumber_vertexesTmp.shape[0], 4), 'i')
        localTestSrcRWGNumber_vertexes[:, 0] = localTestSrcRWGNumber_oppVertexesTmp[:, 0]
        localTestSrcRWGNumber_vertexes[:, 1:3] = localTestSrcRWGNumber_vertexesTmp
        localTestSrcRWGNumber_vertexes[:, 3] = localTestSrcRWGNumber_oppVertexesTmp[:, 1]

        localTestSrcRWGNumber_nodes = localTestSrcRWGNumber_vertexes.flat[:]
        tmpDictionary = {}
        for i in range(len(localTestSrcRWGNumber_nodes)):
            node = localTestSrcRWGNumber_nodes[i]
            if tmpDictionary.has_key(node):
                tmpDictionary[node].append(i)
            else:
                tmpDictionary[node] = [i]
        N_nodes = len(tmpDictionary)

        nodesCoord = zeros((N_nodes, 3), 'd')
        index = 0
        for key, val in tmpDictionary.iteritems():
            nodesCoord[index] = self.vertexes_coord[key]
            for v in val:
                localTestSrcRWGNumber_nodes[v] = index
            index += 1

        cubeNeighborsIndexes = self.cubesNeighborsIndexes[cubeNumber]
        N_neighbors = len(cubeNeighborsIndexes)
        rCubeCenter = self.cubes_centroids[cubeNumber]
        # cubeIntArrays
        cubeIntArrays = zeros(5 + N_RWG_src + 2*N_RWG_src + 4*N_RWG_src + N_RWG_test + N_RWG_src + N_neighbors, 'i')
        cubeIntArrays[:5] = array([N_RWG_test, N_RWG_src, N_neighbors, N_nodes, self.S], 'i')
        
        startIndex = 5
        stopIndex = startIndex + N_RWG_src
        cubeIntArrays[startIndex:stopIndex] = array(list_of_src_RWGsNumbers, 'i')
        
        startIndex = stopIndex
        stopIndex = startIndex + 2*N_RWG_src
        cubeIntArrays[startIndex:stopIndex] = localTestSrcRWGNumber_signedTriangles.flat[:]
        
        startIndex = stopIndex
        stopIndex = startIndex + 4*N_RWG_src
        cubeIntArrays[startIndex:stopIndex] = localTestSrcRWGNumber_nodes

        startIndex = stopIndex
        stopIndex = startIndex + N_RWG_test
        cubeIntArrays[startIndex:stopIndex] = localTestRWGNumber_CFIE_OK
        
        startIndex = stopIndex
        stopIndex = startIndex + N_RWG_src
        cubeIntArrays[startIndex:stopIndex] = localSrcRWGNumber_M_CURRENT_OK
        
        startIndex = stopIndex
        stopIndex = startIndex + N_neighbors
        cubeIntArrays[startIndex:stopIndex] = array(cubeNeighborsIndexes, 'i')
        # cubeDoubleArrays
        cubeDoubleArrays = zeros(nodesCoord.shape[0]*3 + 3, 'd')
        cubeDoubleArrays[:nodesCoord.shape[0]*3] = nodesCoord.flat[:]
        cubeDoubleArrays[nodesCoord.shape[0]*3:] = rCubeCenter
        return cubeIntArrays, cubeDoubleArrays

    def write_cubes(self):
        write_cubes(os.path.join(self.path, self.targetName) + '.cubes', self.cubes_centroids, self.a)


class CubeClass:
    def __init__(self):
        pass
    
    def writeIntDoubleArraysToFile(self, pathToSaveTo, cubeNumber):
        writeBlitzArrayToDisk(self.cubeIntArrays, os.path.join(pathToSaveTo, str(cubeNumber) + "_IntArrays.txt"))
        writeBlitzArrayToDisk(self.cubeDoubleArrays, os.path.join(pathToSaveTo, str(cubeNumber) + "_DoubleArrays.txt"))

    def setIntDoubleArraysFromFile(self, pathToReadFrom, cubeNumber):
        cubeIntArrays = read1DBlitzArrayFromDisk(os.path.join(pathToReadFrom, str(cubeNumber) + "_IntArrays.txt"), 'i')
        self.N_cubeIntArrays = len(cubeIntArrays)
        self.N_RWG_test = cubeIntArrays[0]
        self.N_RWG_src = cubeIntArrays[1]
        self.N_neighbors = cubeIntArrays[2]
        self.N_nodes = cubeIntArrays[3]
        self.S = cubeIntArrays[4]
        # now we copy the arrays
        startIndex = 5
        #stopIndex = startIndex + self.N_RWG_test
        #self.test_RWGsNumbers = cubeIntArrays[startIndex:stopIndex]
    
        #startIndex = stopIndex
        stopIndex = startIndex + self.N_RWG_src
        self.testSrc_RWGsNumbers = cubeIntArrays[startIndex:stopIndex]
    
        startIndex = stopIndex
        stopIndex = startIndex + self.N_RWG_src * 2
        self.localTestSrcRWGNumber_signedTriangles = reshape(cubeIntArrays[startIndex:stopIndex], (-1, 2))
        
        startIndex = stopIndex
        stopIndex = startIndex + self.N_RWG_src * 4
        self.localTestSrcRWGNumber_nodes = reshape(cubeIntArrays[startIndex:stopIndex], (-1, 4))
    
        startIndex = stopIndex
        stopIndex = startIndex + self.N_RWG_test
        self.localTestRWGNumber_CFIE_OK = cubeIntArrays[startIndex:stopIndex]
        
        startIndex = stopIndex
        stopIndex = startIndex + self.N_RWG_src
        self.localSrcRWGNumber_M_CURRENT_OK = cubeIntArrays[startIndex:stopIndex]
        
        startIndex = stopIndex
        stopIndex = startIndex + self.N_neighbors
        self.cubeNeighborsIndexes = cubeIntArrays[startIndex:stopIndex]
        
        # the double arrays
        self.N_cubeDoubleArrays = self.N_nodes * 3 + 3
        #cubeDoubleArrays = zeros(self.N_cubeDoubleArrays, 'd')
        cubeDoubleArrays = read1DBlitzArrayFromDisk(os.path.join(pathToReadFrom, str(cubeNumber) + "_DoubleArrays.txt"), 'd')
        
        startIndex = 0
        stopIndex = startIndex + self.N_nodes * 3
        self.nodesCoord = reshape(cubeDoubleArrays[startIndex:stopIndex], (-1, 3))
        
        startIndex = stopIndex
        stopIndex = startIndex + 3
        self.rCubeCenter = cubeDoubleArrays[startIndex:stopIndex]
   
    def getIntArrays(self):
        return self.N_RWG_test, self.N_RWG_src, self.N_neighbors, self.N_nodes, self.S, self.testSrc_RWGsNumbers, self.localTestSrcRWGNumber_signedTriangles, self.localTestSrcRWGNumber_nodes, self.localTestRWGNumber_CFIE_OK, self.localSrcRWGNumber_M_CURRENT_OK, self.cubeNeighborsIndexes

    def getDoubleArrays(self):
        return self.nodesCoord, self.rCubeCenter


def testScalingMesh(path, targetName, lc, languageForMeshConstruction):
    z_offset = 0.0
    targetDimensions_scaling_factor = 1.0
    NRWG = []
    time_meshing, time_reading, time_edges_classification, time_reordering_normals, time_effective_RWG_functions_computation = [], [], [], [], []
    for i in lc:
        write_geo(path, targetName, 'lc', i)
        t0 = time.time()
        executeGmsh(path, targetName, 0)
        time_meshing.append(time.time()-t0)
        target_mesh = MeshClass(path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction)
        target_mesh.constructFromGmshFile()
        time_reading.append(target_mesh.time_reading)
        time_edges_classification.append(target_mesh.time_edges_classification)
        time_reordering_normals.append(target_mesh.time_reordering_normals)
        time_effective_RWG_functions_computation.append(target_mesh.time_effective_RWG_functions_computation)
        NRWG.append(target_mesh.N_RWG)
    from pylab import plot, rc, subplot, xlabel, ylabel, legend, xticks, yticks, grid, gca, setp, show, title
    #rc('text', usetex=True)
    FontSize=20
    LineWidth=1
    plot(array(NRWG), array(time_meshing), 'ko-', array(NRWG), array(time_reading), 'bo-', array(NRWG), array(time_edges_classification), 'ro-', array(NRWG), array(time_reordering_normals), 'go-', array(NRWG), array(time_effective_RWG_functions_computation), 'yo-')
    legend([r'meshing',r'reading',r'edges classification',r'reordering normals',r'effective RWG functions computation'])
    grid(True)
    show()


if __name__=="__main__":
    path = './geo'
    targetName = 'vf128'
    f = 0.8e9
    write_geo(path, targetName, 'lc', c/f/10.0)
    write_geo(path, targetName, 'lx', 0.1)
    write_geo(path, targetName, 'ly', 0.1)
    write_geo(path, targetName, 'lz', 0.1)
    if 1:
        executeGmsh(path, targetName, 0)
        z_offset = 0.0
        targetDimensions_scaling_factor = 1.0
        languageForMeshConstruction = "C++"
        meshFormat = 'GMSH' 
        meshFileTermination = '.msh'
        target_mesh = MeshClass(path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction, meshFormat, meshFileTermination)
        target_mesh.constructFromGmshFile()
        a = c/f * 0.25
        target_mesh.cubes_data_computation(a)
        target_mesh.write_cubes()
        #target_mesh.saveToDisk("./geo")
        #target_mesh2 = MeshClass(path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction, meshFormat, meshFileTermination)
        #target_mesh2.constructFromSavedArrays("./geo")
        #print "delta_gap =", target_mesh.DELTA_GAP
    if 0:
        #lc = [0.0005]
        lc = [.02, 0.01, 0.0075, 0.005, 0.0025, 0.001, 0.00075]
        testScalingMesh(path, targetName, lc)




