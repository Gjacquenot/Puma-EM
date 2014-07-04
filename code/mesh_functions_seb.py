import os.path, sys, time
from scipy import zeros, ones, arange, array, take, reshape, sort, argsort, put, sum, compress, nonzero, prod, floor, mean, sqrt, dot, arccos
from read_mesh import read_mesh_GMSH_1, read_mesh_GMSH_2
from PyGmsh import executeGmsh, write_geo, findParameter, findParameterValue
from EM_constants import *
import copy

def triangles_centroids_computation(vertexes_coord, triangle_vertexes):
    """guess what? computes the centroids of the triangles"""
    triangles_centroids = take( vertexes_coord, triangle_vertexes[:, 0], axis=0 )
    triangles_centroids += take( vertexes_coord, triangle_vertexes[:, 1], axis=0 )
    triangles_centroids += take( vertexes_coord, triangle_vertexes[:, 2], axis=0 )
    triangles_centroids /= 3.0
    return triangles_centroids

def edges_computation(triangle_vertexes, vertexes_coord):
    """This function builds the edges matrix from the triangles"""

    T = triangle_vertexes.shape[0] # T is the number of triangles
    E = T * 3 # there are 3 edges per triangles
    edges_vertexes = zeros((E, 3), 'i') # edges_vertexes[i, :] -> [v_start, v_end, v_opposite]
    # the edge kind is related to its number of triangles it belongs to.
    # 1 is a physical border, 2 is a normal RWG, and 3 or more is a junction

    t0 = time.clock()
    print("    construction of edges_vertexes...")
    sys.stdout.flush()
    # we first construct a flattened view of edges_vertexes, such that
    # all edges corresponding to 1 triangle are on the same line of the array view
    flat_edges_vertexes = edges_vertexes.reshape((T, -1))

    # we then assign efficiently thanks to the flattened view:
    # edge0 = r0 -> r1, opposite = r2
    # edge1 = r1 -> r2, opposite = r0
    # edge2 = r2 -> r0, opposite = r1
    flat_edges_vertexes[:, 0] = flat_edges_vertexes[:, 5] = flat_edges_vertexes[:, 7] = triangle_vertexes[:, 0]
    flat_edges_vertexes[:, 1] = flat_edges_vertexes[:, 3] = flat_edges_vertexes[:, 8] = triangle_vertexes[:, 1]
    flat_edges_vertexes[:, 2] = flat_edges_vertexes[:, 4] = flat_edges_vertexes[:, 6] = triangle_vertexes[:, 2]

    # For assigning a "number" and a "kind" to an edge (a set of two vertexes),
    # we need to find how many occurrences it has in "edges_vertexes". A costful
    # operation, because for each edge we have to compare the edge to all other
    # that appear in the "edges_vertexes" array.
    #
    # However, if we want efficiency, we can create an array of edges sorted as follows:
    # 1) sort the elements of the first two columns alongside dimension 1; 
    # 2) sort the pairs following their 1st element alongside dimension 0;
    # 3) sort the pairs following their 2nd element alongside dimension 0,
    #    but with keeping the order of their first element.
    # In such an array, all occurrences of an edge would be adjacent to each other. So, for
    # counting the occurrences of an edge, one should only count its similar
    # neighbors, thereby greatly reducing the computational cost of the algorithm

    # we now construct an array of sorted "edges_vertexes" following dimension 1 (the columns) 
    col_sorted_e_v = sort(edges_vertexes[:, :2], 1, kind='mergesort')
    #edges_opp_vertexes = edges_vertexes[:, 2]
    del edges_vertexes
    print("time = " + str(time.clock() - t0))

    t0 = time.clock()
    print("    construction of edgeNumber_triangles...")
    sys.stdout.flush()
    edgeNumber_triangles, edgeNumber_vertexes = compute_edgeNumber_triangles(col_sorted_e_v)
    print("    cumulated time = ", str(time.clock() - t0))
    # construction of triangle_adjacentTriangles matrix
    t0 = time.clock()
    triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = compute_triangle_adjacentTriangles(T, edgeNumber_triangles)
    print("time = " + str(time.clock() - t0))
    return edgeNumber_vertexes.astype('i'), edgeNumber_triangles, triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction

def compute_edgeNumber_triangles(col_sorted_e_v):
    # Once the elements of the first 2 columns have been sorted alongside dimension 1,
    # we construct a 1-D array of real numbers, with:
    # 1) the entire part formed by the numbers of the first column
    # 2) the decimal part formed by the numbers of the second column
    E = col_sorted_e_v.shape[0]
    max_decimal = max(col_sorted_e_v[:, 1]) # the maximum element of "col_sorted_e_v[:, 1]"
    X = 10.0
    while X < max_decimal: # we look for smallest "X" such that "1eX > max_decimal"
        X *= 10.0
    decimal_e_v = col_sorted_e_v[:, 0] + col_sorted_e_v[:, 1]/X
    ind_sorted_e_v = argsort(decimal_e_v, kind='mergesort')
    sorted_decimal_e_v = take(decimal_e_v, ind_sorted_e_v, axis=0)

    diff = ones(sorted_decimal_e_v.shape[0], 'd')
    diff[1:] = abs(sorted_decimal_e_v[1:] - sorted_decimal_e_v[:-1])
    indexesEqualPreceding = compress(diff==0.0,arange(len(diff)),axis=0)
    del diff, decimal_e_v, sorted_decimal_e_v
    t0 = time.clock()
    print("        research of the same edges...")
    sys.stdout.flush()
    indexesEqualEdges = {}
    j = 0
    for i in indexesEqualPreceding:
        indexesEqualEdges[j] = [int(ind_sorted_e_v[i-1]), int(ind_sorted_e_v[i])]
        if j>0:
            if indexesEqualEdges[j][0]==indexesEqualEdges[j-1][-1]:
                # in this case we have a junction: 3 or more triangles that
                # have a common edge. We then merge the two lists and
                # there is no need for increasing the index j
                indexesEqualEdges[j-1] += indexesEqualEdges[j][1:]
                del indexesEqualEdges[j]
            else:
                j += 1
        else:
            j += 1
    print("time = " + str(time.clock() - t0))

    # edge numbering
    t0 = time.clock()
    print("        numbering of the edges...")
    sys.stdout.flush()
    edgeNumber_vertexes = ones((len(indexesEqualEdges), 2), 'i') * -1
    number = 0
    for key, value in indexesEqualEdges.items():
        # All occurrences of an "edge" receive the same "edge_number".
        #for i in value:
        #    edges_new_numbers[i] = number
        edgeNumber_vertexes[number] = col_sorted_e_v[value[0]]
        number += 1
    max_edges_numbers = number
    print("time = " + str(time.clock() - t0))

    # construction of edgeNumber_triangles
    t0 = time.clock()
    print("        construction of edgeNumber_triangles...")
    sys.stdout.flush()
    edgeNumber_triangles = indexesEqualEdges
    for key, value in edgeNumber_triangles.items():
        value_mod_3 = [int(x/3) for x in value]
        edgeNumber_triangles[key] = value_mod_3
    print("time = " + str(time.clock() - t0))
    return edgeNumber_triangles, edgeNumber_vertexes

def compute_triangle_adjacentTriangles(T, edgeNumber_triangles):
    print("    construction of triangle_adjacentTriangles...")
    sys.stdout.flush()
    triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = {}, {}
    # initialisation of the dictionaries
    for i in range(T):
        triangle_adjacentTriangles[i] = []
    # filling of the dictionaries
    for adjacent_triangles in edgeNumber_triangles.values():
        N_adj_triangles = len(adjacent_triangles)
        IS_JUNCTION = (N_adj_triangles>2)
        for tn in adjacent_triangles:
            listToAdd = [tj for tj in adjacent_triangles if tj!=tn]
            triangle_adjacentTriangles[tn] += listToAdd
            if IS_JUNCTION:
                if tn in is_triangle_adjacentTriangles_via_junction:
                    is_triangle_adjacentTriangles_via_junction[tn] += listToAdd #[IS_JUNCTION] * (N_adj_triangles-1)
                else:
                    is_triangle_adjacentTriangles_via_junction[tn] = listToAdd
    return triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction

def RWGNumber_signedTriangles_computation(edgeNumber_triangles, edgeNumber_vertexes, triangles_surfaces, is_closed_surface, triangle_vertexes, vertexes_coord):
    # we now want to "intelligently" get rid of the junctions, that means, create new RWGs when needed.
    # The following is correct for metal-metal junctions
    print("    computation of RWG to triangles relations...")
    sys.stdout.flush()
    t5 = time.clock()
    N_edges = len(edgeNumber_triangles)
    # RWGNumber_signedTrianglesTmp_1 is an array of fixed size N_edges. It will be equal to
    # edgeNumber_triangles if there are no junctions. If there are junctions, RWGNumber_signedTrianglesTmp_2 
    # will be non-empty. This is for limiting the memory requirements on this part of the code.
    RWGNumber_signedTrianglesTmp_1, RWGNumber_signedTrianglesTmp_2 = zeros((N_edges, 2), 'i'), {}
    RWGNumber_edgeNumber = range(N_edges)
    for index in range(N_edges):
        numberOfTriangles = len(edgeNumber_triangles[index])
        if numberOfTriangles==2:
            RWGNumber_signedTrianglesTmp_1[index] = array(edgeNumber_triangles[index], 'i')
        else: # then we have a junction
            # we compute the edge unit vector: zHatEdge
            n0, n1 = edgeNumber_vertexes[index, 0], edgeNumber_vertexes[index, 1]
            r0, r1 = vertexes_coord[n0], vertexes_coord[n1]
            r1_r0 = r1 - r0
            zHatEdge = r1_r0 / sqrt(dot(r1_r0, r1_r0))
            xHatEdge, yHatEdge = zeros(3, 'd'), zeros(3, 'd')
            triangles = edgeNumber_triangles[index]
            triangles_angles = zeros(numberOfTriangles, 'd')
            # construction of the edge local coordinate system: it is based on the first triangle
            # this is for sorting the triangles wrt their respective angles
            # because we cannot have a RWG which is bisected by another RWG
            for i in range(numberOfTriangles):
                tr = triangles[i]
                tr_nodes = triangle_vertexes[tr]
                for n2 in tr_nodes:
                    if (n2!= n0) and (n2!= n1):
                        break
                r2 = vertexes_coord[n2]
                r2_r0 = r2 - r0
                r2_r0_perpendicular = r2_r0 - dot(r2_r0, zHatEdge) * zHatEdge
                r2_r0_perpendicularHat = r2_r0_perpendicular / sqrt(dot(r2_r0_perpendicular, r2_r0_perpendicular))
                if i==0:
                    xHatEdge = r2_r0_perpendicularHat[:] # the first triangle will be the reference: will define the xHat vector
                    yHatEdge[0] = zHatEdge[1] * xHatEdge[2] - zHatEdge[2] * xHatEdge[1]
                    yHatEdge[1] = zHatEdge[2] * xHatEdge[0] - zHatEdge[0] * xHatEdge[2]
                    yHatEdge[2] = zHatEdge[0] * xHatEdge[1] - zHatEdge[1] * xHatEdge[0]
                    triangles_angles[0] = 0.0
                else:
                    Cos = dot(r2_r0_perpendicularHat, xHatEdge)
                    Sin = dot(r2_r0_perpendicularHat, yHatEdge)
                    if Sin>=0:
                        triangles_angles[i] = arccos(Cos)
                    else:
                        triangles_angles[i] = 2*pi - arccos(Cos)
            # we now sort the triangles wrt their respective position wrt xHatEdge
            ind_sortedTriangles = argsort(triangles_angles, kind='mergesort') 
            sortedTriangles = take(triangles, ind_sortedTriangles)
            sortedTriangleSurfaces = triangles_surfaces[sortedTriangles]
            # we now form all the possible RWGs with the sorted triangles. Normally none of these RWGs can be bisected by a triangle now
            possibleTrianglesPairsForRWGs = []
            for i in range(1, numberOfTriangles):
                if (sortedTriangleSurfaces[i] != sortedTriangleSurfaces[i-1]):
                    possibleTrianglesPairsForRWGs.append([sortedTriangles[i-1], sortedTriangles[i]])
            # the last possibility of RWG: between the last and first triangle
            if len(possibleTrianglesPairsForRWGs)<numberOfTriangles-1: # if our list is too small
                if (sortedTriangleSurfaces[numberOfTriangles-1] != sortedTriangleSurfaces[0]): # if the triangles belong to a different surface
                    possibleTrianglesPairsForRWGs.append([sortedTriangles[numberOfTriangles-1], sortedTriangles[0]])
            for i in range(1, len(possibleTrianglesPairsForRWGs)): # if more than one, we have a junction
                key = len(RWGNumber_signedTrianglesTmp_2)
                RWGNumber_signedTrianglesTmp_2[key] = possibleTrianglesPairsForRWGs[i]
                RWGNumber_edgeNumber.append(index)
            # we don't need the whole list of edges for current line of RWGNumber_signedTriangles
            RWGNumber_signedTrianglesTmp_1[index] = array(possibleTrianglesPairsForRWGs[0], 'i')
    # we can finally construct the final RWGNumber_signedTriangles array
    N_RWG = RWGNumber_signedTrianglesTmp_1.shape[0] + len(RWGNumber_signedTrianglesTmp_2)
    RWGNumber_signedTriangles = zeros((N_RWG, 2), 'i')
    RWGNumber_signedTriangles[:N_edges, :] = RWGNumber_signedTrianglesTmp_1
    index = N_edges
    for key, value in RWGNumber_signedTrianglesTmp_2.items():
        RWGNumber_signedTriangles[index] = array(RWGNumber_signedTrianglesTmp_2[key], 'i')
        index += 1
    if not (index==N_RWG):
        print("Error at the end of RWGNumber_signedTriangles_computation. Exiting")
        sys.exit(1)
    # computation of RWGNumber_edgeVertexes
    #RWGNumber_edgeVertexes = (take(edgeNumber_vertexes, array(RWGNumber_edgeNumber, 'i'), axis=0)).astype('i')
    RWGNumber_edgeVertexes = zeros((N_RWG, 2), 'i')
    for i in range(N_RWG):
        t0, t1 = RWGNumber_signedTriangles[i, 0], RWGNumber_signedTriangles[i, 1]
        s0, s1 = triangles_surfaces[t0], triangles_surfaces[t1]
        t = t0
        if not s0==s1:
            if is_closed_surface[s0]==1:
                t = t0
            elif is_closed_surface[s1]==1:
                t = t1
            else:
                t = t0
        e0, e1 = edgeNumber_vertexes[RWGNumber_edgeNumber[i], 0], edgeNumber_vertexes[RWGNumber_edgeNumber[i], 1]
        n0, n1, n2 = triangle_vertexes[t, 0], triangle_vertexes[t, 1], triangle_vertexes[t, 2]
        if ((e0==n0) and (e1==n1)) or ((e0==n1) and (e1==n2)) or ((e0==n2) and (e1==n0)):
            RWGNumber_edgeVertexes[i, 0] = e0
            RWGNumber_edgeVertexes[i, 1] = e1
        else:
            RWGNumber_edgeVertexes[i, 0] = e1
            RWGNumber_edgeVertexes[i, 1] = e0
    print("   time = " + str(time.clock() - t5))
    return RWGNumber_signedTriangles.astype('i'), RWGNumber_edgeVertexes.astype('i'), N_edges, N_RWG

def RWGNumber_oppVertexes_computation(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, triangle_vertexes):
    print("    computation of RWG opposite vertexes...")
    sys.stdout.flush()
    t5 = time.clock()
    N_RWG = RWGNumber_signedTriangles.shape[0]
    RWGNumber_oppVertexes = zeros((N_RWG, 2), 'i')
    for j in range(N_RWG):
        t0 = RWGNumber_signedTriangles[j, 0]
        t1 = abs(RWGNumber_signedTriangles[j, 1])
        edgeVertexes = RWGNumber_edgeVertexes[j]
        t0Vertexes = triangle_vertexes[t0]
        t1Vertexes = triangle_vertexes[t1]
        for i in t0Vertexes:
            if (i!=edgeVertexes[0]) and (i!=edgeVertexes[1]):
                break
        RWGNumber_oppVertexes[j, 0] = i
        for i in t1Vertexes:
            if (i!=edgeVertexes[0]) and (i!=edgeVertexes[1]):
                break
        RWGNumber_oppVertexes[j, 1] = i
    print("   time = " + str(time.clock() - t5))
    return RWGNumber_oppVertexes.astype('i')

def compute_RWGNumber_edgeCentroidCoord(vertexes_coord, RWGNumber_edgeVertexes):
    edgeCentroidCoord = take( vertexes_coord, RWGNumber_edgeVertexes[:, 0], axis=0 )
    edgeCentroidCoord += take( vertexes_coord, RWGNumber_edgeVertexes[:, 1], axis=0 )
    return edgeCentroidCoord/2.0

def compute_RWGNumber_edgeLength(vertexes_coord, RWGNumber_edgeVertexes):
    RWGNumber_length = take( vertexes_coord, RWGNumber_edgeVertexes[:, 0], axis=0 )
    RWGNumber_length -= take( vertexes_coord, RWGNumber_edgeVertexes[:, 1], axis=0 )
    RWGNumber_length = sqrt(sum((RWGNumber_length)**2, axis=1))
    return RWGNumber_length.astype('d')

def compute_RWG_meanEdgeLength(vertexes_coord, RWGNumber_edgeVertexes, stride):
    r0 = take( vertexes_coord, RWGNumber_edgeVertexes[::stride, 0], axis=0 )
    r1 = take( vertexes_coord, RWGNumber_edgeVertexes[::stride, 1], axis=0 )
    RWGNumber_length = sqrt(sum((r0 - r1)**2, axis=1))
    return sum(RWGNumber_length)/RWGNumber_length.shape[0]


def edgeNumber_triangles_indexes(list_of_edges_numbers, RWGNumber_signedTriangles):
    """This function returns a 1-D array of the indexes of the triangles corresponding
    to a 1-D array of edges_numbers. This function is important for creating lists of triangles
    that will participate to the MoM, given a particular criterium concerning the edges."""
    indexes_of_triangles_tmp1 = take(RWGNumber_signedTriangles, list_of_edges_numbers, axis=0).flat
    indexes_of_triangles_tmp2 = sort(indexes_of_triangles_tmp1, kind='mergesort')
    indexes_of_triangles_to_take = ones(indexes_of_triangles_tmp2.shape[0], 'i')
    indexes_of_triangles_to_take[1:] = indexes_of_triangles_tmp2[1:] - indexes_of_triangles_tmp2[:-1]
    indexes_of_triangles = compress(indexes_of_triangles_to_take != 0, indexes_of_triangles_tmp2)
    return indexes_of_triangles.astype('i')

def change_triangle_circulation(t0, t1, triangle_vertexes):
    """This function reorders the nodes of triangle t1 such that the common edge
    between reference triangle t0 and t1 is parcouru following one direction and
    then following the opposite. In this way we will have compatible normals"""
    # nodes of reference triangle
    node00 = triangle_vertexes[t0, 0]
    node01 = triangle_vertexes[t0, 1]
    node02 = triangle_vertexes[t0, 2]
    # nodes of triangle to change
    node10 = triangle_vertexes[t1, 0]
    node11 = triangle_vertexes[t1, 1]
    node12 = triangle_vertexes[t1, 2]
    # circulation on these triangles is a succession of paths between nodes: 0 -> 1, 1 -> 2, 2 -> 0
    t0_circ = [[node00, node01], [node01, node02], [node02, node00]]
    t1_circ = [[node10, node11], [node11, node12], [node12, node10]]
    # "coincide_circ": a variable that will tell if the edge is parcouru alike by each triangle.
    # In this case, one will have to reorder the nodes of t1.
    coincide_circ = 0;
    for k in range(3):
        if coincide_circ == 1:
            break
        for l in range(3):
            if t0_circ[k] == t1_circ[l]:
                coincide_circ = 1
                break
    if coincide_circ == 1:
        triangle_vertexes[t1, :] = array([node10, node12, node11])


def reorder_triangle_vertexes(triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction, triangle_vertexes, vertexes_coord):
    """This function is necessary to orientate all the normals of the triangles
    coherently for a given surface. For this purpose it gives a good order of
    appearance of the triangles vertexes. This function acts directly on
    the array 'triangle_vertexes'.
    
    This function also returns triangles_surfaces"""

    print("      reordering triangles for normals coherency...")
    sys.stdout.flush()
    T = len(triangle_adjacentTriangles)
    is_triangle_reordered = [0] * T # create a list of length "T" that tells if the triangle has been reordered
    triangles_surfaces = array([-1] * T, 'i')
    surf_number = -1
    while 0 in is_triangle_reordered: # as long as there is a triangle that has not been reordered
        surf_number += 1
        t_start = is_triangle_reordered.index(0) # the triangle from which we start the growing
        is_triangle_reordered[t_start] = 1
        triangles_surfaces[t_start] = surf_number
        #index = 0
        list_t_to_reorder = []
        for tn in triangle_adjacentTriangles[t_start]:
            is_triangle_not_adjacent_via_junction = True
            if t_start in is_triangle_adjacentTriangles_via_junction:
                if tn in is_triangle_adjacentTriangles_via_junction[t_start]:
                    is_triangle_not_adjacent_via_junction = False
            if (is_triangle_reordered[int(tn)]==0) and (is_triangle_not_adjacent_via_junction):
                list_t_to_reorder.append(tn)
            #index += 1
        # list_t_to_reorder = [tn for tn in triangle_adjacentTriangles[t_start] if is_triangle_reordered[int(tn)]==0] # initialization of the list to process
        list_calling_t = [t_start] * len(list_t_to_reorder)
        while list_t_to_reorder:
            t = list_t_to_reorder.pop()
            calling_t = list_calling_t.pop()
            # we change circulation of t according to calling_t
            change_triangle_circulation(calling_t, t, triangle_vertexes) 
            is_triangle_reordered[int(t)] = 1
            triangles_surfaces[int(t)] = surf_number
            #index = 0
            t_adjacent_triangles = []
            for tn in triangle_adjacentTriangles[t]:
                is_triangle_not_adjacent_via_junction = True
                if t in is_triangle_adjacentTriangles_via_junction:
                    if tn in is_triangle_adjacentTriangles_via_junction[t]:
                        is_triangle_not_adjacent_via_junction = False
                if (is_triangle_reordered[int(tn)]==0) and (is_triangle_not_adjacent_via_junction):
                    t_adjacent_triangles.append(tn)
                #index += 1
            # we extend the list to reorder with the triangles adjacent to t
            list_t_to_reorder.extend(t_adjacent_triangles)
            # we extend the list of "calling" triangles with t
            list_calling_t.extend([int(t)] * len(t_adjacent_triangles))  
    # we loop on the surfaces, because all normals are coherent but maybe not directed outwards closed surfaces...
    print("      redirecting the normals outward...")
    sys.stdout.flush()
    S = surf_number
    triangles_centroids_z = triangles_centroids_computation(vertexes_coord, triangle_vertexes)[:,2]
    for s in range(S+1):
        ind_t_on_s = nonzero(triangles_surfaces == s)[0]
        max_height_centroids_s = max(take(triangles_centroids_z, ind_t_on_s, axis=0))
        ind_max_height_centroids_s_tmp = take(triangles_centroids_z, ind_t_on_s, axis=0).tolist().index(max_height_centroids_s)
        ind_max_height_centroids_s = ind_t_on_s[ind_max_height_centroids_s_tmp]
        r0 = vertexes_coord[triangle_vertexes[ind_max_height_centroids_s, 0]]
        r1 = vertexes_coord[triangle_vertexes[ind_max_height_centroids_s, 1]]
        r2 = vertexes_coord[triangle_vertexes[ind_max_height_centroids_s, 2]]
        r1_r0 = r1 - r0
        r2_r0 = r2 - r0
        triangle_normal = zeros(3, 'd')
        triangle_normal[0] = r1_r0[1] * r2_r0[2] - r1_r0[2] * r2_r0[1]
        triangle_normal[1] = r1_r0[2] * r2_r0[0] - r1_r0[0] * r2_r0[2]
        triangle_normal[2] = r1_r0[0] * r2_r0[1] - r1_r0[1] * r2_r0[0]
        # we now test the normal associated to ind_max_height_centroids_s
        z_hat = array([0, 0, 1], 'd')
        if sum(triangle_normal * z_hat) < 0:
            # we also swap the columns of triangle_vertexes
            for index in ind_t_on_s: # this part does not work with put...
                vertexes = copy.copy(triangle_vertexes[index, :])
                triangle_vertexes[index, 1] = vertexes[2]
                triangle_vertexes[index, 2] = vertexes[1]
    return triangles_surfaces.astype('i')

def is_surface_closed(triangles_surfaces, edgeNumber_triangles):
    """this function determines if a surface is open or closed.
       it also provides relationships between surfaces (linked or not)"""
    S = max(triangles_surfaces)+1
    NUMBER_TRIANGLES_IN_SURFACE = zeros(S, 'i')
    for s in triangles_surfaces:
        NUMBER_TRIANGLES_IN_SURFACE[s] += 1
    connected_surfaces = {}
    # we now count the number of INNER edges for each surface.
    # the edges that are junctions receive a special treatment:
    # only if the edge has two triangles on the given surface, 
    # can it be counted as an inner edge, which will then be 
    # counted in NUMBER_EDGES_IN_SURFACE
    NUMBER_EDGES_IN_SURFACE = zeros(S, 'i')
    for edge_number, triangles_tmp in edgeNumber_triangles.items():
        surfaces_appeared_already = zeros(S, 'i')
        if len(triangles_tmp)>2: # we have a junction here
            for t in triangles_tmp:
                surface = triangles_surfaces[t]
                if surfaces_appeared_already[surface]==0:
                    surfaces_appeared_already[surface] = 1
                else:
                    NUMBER_EDGES_IN_SURFACE[surface] += 1
            surfaces_present = compress(surfaces_appeared_already>0, arange(S))
            if len(surfaces_present)==2:
                s0, s1 = min(surfaces_present), max(surfaces_present)
                if (s0, s1) in connected_surfaces:
                    connected_surfaces[(s0, s1)].append(edge_number)
                else:
                    connected_surfaces[(s1, s0)] = [edge_number]
            else:
                for index1 in arange(len(surfaces_present)):
                    for index2 in arange(index1+1, len(surfaces_present)):
                        s1 = min(surfaces_present[index1], surfaces_present[index2])
                        s2 = max(surfaces_present[index1], surfaces_present[index2])
                        if (s1, s2) in connected_surfaces:
                            connected_surfaces[(s1, s2)].append(edge_number)
                        else:
                            connected_surfaces[(s1, s2)] = [edge_number]
        else:
            surface = triangles_surfaces[triangles_tmp[0]]
            NUMBER_EDGES_IN_SURFACE[surface] += 1

    is_closed_surface = ( (NUMBER_EDGES_IN_SURFACE*2) == (NUMBER_TRIANGLES_IN_SURFACE*3) )
    # we now check for potential closed surfaces: surfaces which can be closed
    # and on which we can therefore apply the CFIE
    potential_closed_surfaces = {}
    for key, item in connected_surfaces.items():
        s0, s1 = key[0], key[1]
        numberEdges0, numberEdges1 = NUMBER_EDGES_IN_SURFACE[s0], NUMBER_EDGES_IN_SURFACE[s1]
        numberTriangles0, numberTriangles1 = NUMBER_TRIANGLES_IN_SURFACE[s0], NUMBER_TRIANGLES_IN_SURFACE[s1]
        if ( numberEdges0 + numberEdges1 + len(item) )*2 == 3*(numberTriangles0 + numberTriangles1):
            potential_closed_surfaces[key] = item
    return is_closed_surface * 1, connected_surfaces, potential_closed_surfaces

def triangles_unnormalized_normals_computation(vertexes_coord, triangle_vertexes, t):
    """This function returns the non-normalized normals of each triangle.
    t is an array of the t indexes to be considered"""
    triangles_normals = zeros((t.shape[0], 3), 'd')
    stride = 10000
    startIndex, stopIndex = 0, min(stride, t.shape[0])
    # it is coded this way for memory optimization: this function is really a memory hog!
    while startIndex<triangles_normals.shape[0]:
        indexes = t[startIndex:stopIndex]
        v0 = take(triangle_vertexes[:, 0], indexes, axis=0)
        v1 = take(triangle_vertexes[:, 1], indexes, axis=0)
        v2 = take(triangle_vertexes[:, 2], indexes, axis=0)
        r0 = take(vertexes_coord, v0, axis=0) # first vertexes of all triangles
        r1_r0 = take(vertexes_coord, v1, axis=0) - r0
        r2_r0 = take(vertexes_coord, v2, axis=0) - r0
        triangles_normals[indexes, 0] = r1_r0[:, 1] * r2_r0[:, 2] - r1_r0[:, 2] * r2_r0[:, 1]
        triangles_normals[indexes, 1] = r1_r0[:, 2] * r2_r0[:, 0] - r1_r0[:, 0] * r2_r0[:, 2]
        triangles_normals[indexes, 2] = r1_r0[:, 0] * r2_r0[:, 1] - r1_r0[:, 1] * r2_r0[:, 0]
        startIndex = stopIndex
        stopIndex = min(stopIndex + stride, t.shape[0])
    return triangles_normals.astype('d')


def triangles_areas_normals_computation(vertexes_coord, triangle_vertexes, triangles_surfaces):
    """this function"""
    T = triangle_vertexes.shape[0]
    S = max(triangles_surfaces)
    triangles_normals = triangles_unnormalized_normals_computation(vertexes_coord, triangle_vertexes, arange(T).astype('i')).astype('d')
    norm_triangles_normals = reshape(sqrt(sum(triangles_normals**2, 1)), (T,1))
    triangles_areas = norm_triangles_normals/2.0
    triangles_normals /= norm_triangles_normals
    return triangles_areas.astype('d'), triangles_normals.astype('d')


def write_normals(name, triangles_centroids, triangles_normals, triangles_surfaces, surface):
    """function that writes the normals of a given surface to a file readable by GMSH for viewing."""
    f = open(name, 'w')
    f.write('View "normals of surfaces xxx" {\n')
    T = triangles_surfaces.shape[0]
    for k in range(T):
        write_condition = (triangles_surfaces[k] == surface) or (surface == -1) # if surface == -1, we write all the normals
        if write_condition:
            string_to_write = 'VP(' + str(triangles_centroids[k, :].tolist())[1:-1] + ')'
            string_to_write += '{' + str(triangles_normals[k,:].tolist())[1:-1] + '};\n'
            f.write(string_to_write)
    f.write('};\n')
    f.close()


def divide_triangles(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, reordered_triangle_vertexes, vertexes_coord):
    N_RWG = RWGNumber_signedTriangles.shape[0]
    V = vertexes_coord.shape[0]
    T = reordered_triangle_vertexes.shape[0]
    divided_triangles_vertexes = zeros((T,7),'int32')
    # in divided_triangles_vertexes, the column indexes:
    # 0, 2, 4 correspond to the nodes 0, 1, 2 of the triangles_vertexes
    # 1, 3, 5 correspond to the midpoints of the edges n01, n12, n20 of the triangles
    # 6 corresponds to the r_grav of the triangle 
    divided_triangles_vertexes[:, 0] = reordered_triangle_vertexes[:, 0]
    divided_triangles_vertexes[:, 2] = reordered_triangle_vertexes[:, 1]
    divided_triangles_vertexes[:, 4] = reordered_triangle_vertexes[:, 2]
    divided_triangles_vertexes[:, 6] = arange(V,V+T)
    # first we take care of the points located on RWG edges
    for i in range(N_RWG):
        e0 = RWGNumber_edgeVertexes[i, 0]
        e1 = RWGNumber_edgeVertexes[i, 1]
        edge = [e0, e1]
        edge.sort()
        number_of_edge_centroid = V + T + i
        for j in range(2):
            t = RWGNumber_signedTriangles[i, j]
            n0 = reordered_triangle_vertexes[t, 0]
            n1 = reordered_triangle_vertexes[t, 1]
            n2 = reordered_triangle_vertexes[t, 2]
            edges_triangle = [[n0,n1], [n1,n2], [n2,n0]]
            for e in edges_triangle:
                e.sort()
            index = edges_triangle.index(edge)
            if divided_triangles_vertexes[t, index*2+1]==0:
                divided_triangles_vertexes[t, index*2+1] = number_of_edge_centroid

    # then we take care of the points located on triangle edges that are not RWGs
    number_of_edge_centroid = V + T + N_RWG
    for t in range(T):
        for j in range(3):
            if (divided_triangles_vertexes[t, j*2+1] == 0):
                divided_triangles_vertexes[t, j*2+1] = number_of_edge_centroid
                number_of_edge_centroid += 1

    # that's all for the barycentric division of the triangles.
    return divided_triangles_vertexes, number_of_edge_centroid

def create_barycentric_triangles(divided_triangles_vertexes, vertexes_coord, MAX_V):
    T = divided_triangles_vertexes.shape[0]
    V = vertexes_coord.shape[0]
    vertexes_coord_barycentric = zeros((MAX_V, 3),'d')
    vertexes_coord_barycentric[:V,:] = vertexes_coord
    barycentric_triangles_vertexes = zeros((T*6,3), 'int32')
    for t in range(T):
        # REMINDER: in divided_triangles_vertexes, the column indexes:
        # 0, 2, 4 correspond to the nodes 0, 1, 2 of the triangles_vertexes
        # 1, 3, 5 correspond to the midpoints of the edges n01, n12, n20 of the triangles
        # 6 corresponds to the r_grav of the triangle
        # draw a triangle and bary-divide it and number points with above indexes to derive the following 
        n0 = divided_triangles_vertexes[t, 0]
        n1 = divided_triangles_vertexes[t, 1]
        n2 = divided_triangles_vertexes[t, 2]
        n3 = divided_triangles_vertexes[t, 3]
        n4 = divided_triangles_vertexes[t, 4]
        n5 = divided_triangles_vertexes[t, 5]
        n6 = divided_triangles_vertexes[t, 6]

        barycentric_triangles_vertexes[t*6+0, :] = [n0, n1, n6]
        barycentric_triangles_vertexes[t*6+1, :] = [n2, n6, n1]
        barycentric_triangles_vertexes[t*6+2, :] = [n2, n3, n6]
        barycentric_triangles_vertexes[t*6+3, :] = [n4, n6, n3]
        barycentric_triangles_vertexes[t*6+4, :] = [n4, n5, n6]
        barycentric_triangles_vertexes[t*6+5, :] = [n0, n6, n5]

        vertexes_coord_barycentric[n1,:] = (vertexes_coord[n0,:] + vertexes_coord[n2,:])/2.0
        vertexes_coord_barycentric[n3,:] = (vertexes_coord[n2,:] + vertexes_coord[n4,:])/2.0
        vertexes_coord_barycentric[n5,:] = (vertexes_coord[n4,:] + vertexes_coord[n0,:])/2.0
        vertexes_coord_barycentric[n6,:] = (vertexes_coord[n0,:] + vertexes_coord[n2,:] + vertexes_coord[n4,:])/3.0

    return vertexes_coord_barycentric, barycentric_triangles_vertexes

def create_barycentric_RWGs(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, divided_triangles_vertexes, barycentric_triangles_vertexes):
    T = divided_triangles_vertexes.shape[0]
    T_bary = 6*T
    N_RWG = RWGNumber_signedTriangles.shape[0]
    barycentric_RWGNumber_signedTriangles = zeros((T_bary+2*N_RWG, 2),'int32')
    barycentric_RWGNumber_edgeVertexes = zeros((T_bary+2*N_RWG, 2),'int32')
    barycentric_RWGNumber_oppVertexes = zeros((T_bary+2*N_RWG, 2),'int32')    
    # first we take care of the 6 barycentric RWGs that are enclosed in each original triangle 
    for t in range(T):
        barycentric_RWGNumber_signedTriangles[6*t + 0, :] = [6*t+0, 6*t+1]
        barycentric_RWGNumber_signedTriangles[6*t + 1, :] = [6*t+1, 6*t+2]
        barycentric_RWGNumber_signedTriangles[6*t + 2, :] = [6*t+2, 6*t+3]
        barycentric_RWGNumber_signedTriangles[6*t + 3, :] = [6*t+3, 6*t+4]
        barycentric_RWGNumber_signedTriangles[6*t + 4, :] = [6*t+4, 6*t+5]
        barycentric_RWGNumber_signedTriangles[6*t + 5, :] = [6*t+5, 6*t+0]

        # REMINDER: in divided_triangles_vertexes, the column indexes:
        # 0, 2, 4 correspond to the nodes 0, 1, 2 of the triangles_vertexes
        # 1, 3, 5 correspond to the midpoints of the edges n01, n12, n20 of the triangles
        # 6 corresponds to the r_grav of the triangle
        # draw a triangle and bary-divide it and number points with above indexes to derive the following 
        barycentric_RWGNumber_edgeVertexes[6*t + 0, :] = [divided_triangles_vertexes[t,1], divided_triangles_vertexes[t,6]]
        barycentric_RWGNumber_edgeVertexes[6*t + 1, :] = [divided_triangles_vertexes[t,2], divided_triangles_vertexes[t,6]]
        barycentric_RWGNumber_edgeVertexes[6*t + 2, :] = [divided_triangles_vertexes[t,3], divided_triangles_vertexes[t,6]]
        barycentric_RWGNumber_edgeVertexes[6*t + 3, :] = [divided_triangles_vertexes[t,4], divided_triangles_vertexes[t,6]]
        barycentric_RWGNumber_edgeVertexes[6*t + 4, :] = [divided_triangles_vertexes[t,5], divided_triangles_vertexes[t,6]]
        barycentric_RWGNumber_edgeVertexes[6*t + 5, :] = [divided_triangles_vertexes[t,0], divided_triangles_vertexes[t,6]]

        barycentric_RWGNumber_oppVertexes[6*t + 0, :] = [divided_triangles_vertexes[t,0], divided_triangles_vertexes[t,2]]
        barycentric_RWGNumber_oppVertexes[6*t + 1, :] = [divided_triangles_vertexes[t,1], divided_triangles_vertexes[t,3]]
        barycentric_RWGNumber_oppVertexes[6*t + 2, :] = [divided_triangles_vertexes[t,2], divided_triangles_vertexes[t,4]]
        barycentric_RWGNumber_oppVertexes[6*t + 3, :] = [divided_triangles_vertexes[t,3], divided_triangles_vertexes[t,5]]
        barycentric_RWGNumber_oppVertexes[6*t + 4, :] = [divided_triangles_vertexes[t,4], divided_triangles_vertexes[t,0]]
        barycentric_RWGNumber_oppVertexes[6*t + 5, :] = [divided_triangles_vertexes[t,5], divided_triangles_vertexes[t,1]]

    # we then create the 2 barycentric RWGs that are associated with each original RWG
    for i in range(N_RWG):
        e0 = RWGNumber_edgeVertexes[i, 0]
        e1 = RWGNumber_edgeVertexes[i, 1]

        t0 = RWGNumber_signedTriangles[i, 0]
        t1 = RWGNumber_signedTriangles[i, 1]
        nodes_t0 = divided_triangles_vertexes[t0,:].tolist()
        nodes_t1 = divided_triangles_vertexes[t1,:].tolist()
        index_start_node = nodes_t0.index(e0)
        index_end_node = nodes_t0.index(e1)
        index_mid_node = (index_start_node + index_end_node)/2
        if ((index_start_node==0) and (index_end_node==4)) or ((index_start_node==4) and (index_end_node==0)):
            index_mid_node = 5
        # the points forming the barycentric RWGs
        barycentric_RWGNumber_edgeVertexes[T_bary + 2*i, :] = [nodes_t0[index_start_node], nodes_t0[index_mid_node]]
        barycentric_RWGNumber_edgeVertexes[T_bary + 2*i + 1, :] = [nodes_t0[index_mid_node], nodes_t0[index_end_node]]
        # the points opposite the barycentric RWGs
        barycentric_RWGNumber_oppVertexes[T_bary + 2*i, :] = [nodes_t0[6], nodes_t1[6]]
        barycentric_RWGNumber_oppVertexes[T_bary + 2*i + 1, :] = [nodes_t0[6], nodes_t1[6]]
        # the signed barycentric triangles forming the RWG: first the two triangles in t0
        start_node = nodes_t0[index_start_node]
        mid_node = nodes_t0[index_mid_node]
        end_node = nodes_t0[index_end_node]
        for j in range(6):
            t0_j = barycentric_triangles_vertexes[6*t0+j, :].tolist()
            t1_j = barycentric_triangles_vertexes[6*t1+j, :].tolist()
            if (start_node in t0_j) and (mid_node in t0_j):
                barycentric_RWGNumber_signedTriangles[T_bary + 2*i, 0] = 6*t0 + j
            if (mid_node in t0_j) and (end_node in t0_j):
                barycentric_RWGNumber_signedTriangles[T_bary + 2*i + 1, 0] = 6*t0 + j
            if (start_node in t1_j) and (mid_node in t1_j):
                barycentric_RWGNumber_signedTriangles[T_bary + 2*i, 1] = 6*t1 + j
            if (mid_node in t1_j) and (end_node in t1_j):
                barycentric_RWGNumber_signedTriangles[T_bary + 2*i + 1, 1] = 6*t1 + j

    return barycentric_RWGNumber_signedTriangles, barycentric_RWGNumber_edgeVertexes, barycentric_RWGNumber_oppVertexes


def create_RWG_to_barycentricRWG(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, divided_triangles_vertexes, vertexes_coord_barycentric):
    # for a given, original RWG, we have 14 participating barycentric RWGs: 6 per triangle plus 2 on the edge.
    # of these 14 barycentric RWGs, only 10 actually participate for rebuilding the original RWG.
    # The 14 can be related to the original RWG as follows:
    # Original RWG i -> 2 original triangles: RWGNumber_signedTriangles[i, :]
    # Original triangle t -> barycentric RWG j: barycentric_RWGNumber_signedTriangles[6*t + j, :]
    # We have then 12 RWGs, but we still lack the two barycentric RWGs that are on the original RWG edge.
    # Original RWG i -> barycentric_RWGNumber_signedTriangles[6*t + 2*i, :] and [6*t + 2*i + 1, :]
    T = divided_triangles_vertexes.shape[0]
    N_RWG = RWGNumber_signedTriangles.shape[0]
    RWG_to_barycentricRWG = zeros((N_RWG, 14), 'int32')
    for i in range(N_RWG):
        t0 = RWGNumber_signedTriangles[i,0]
        t1 = RWGNumber_signedTriangles[i,1]
        for j in range(6):
            RWG_to_barycentricRWG[i, j] = 6*t0 + j
            RWG_to_barycentricRWG[i, j+8] = 6*t1 + j
        RWG_to_barycentricRWG[i, 6] = T*6 + 2*i
        RWG_to_barycentricRWG[i, 7] = T*6 + 2*i + 1

    RWG_to_barycentricRWG_coefficients = zeros((N_RWG, 14), 'd')
    # we already know that the two barycentric RWGs that are on the
    # edge of the RWG have a coefficient of 1 (Andriulli TAP 2008).
    RWG_to_barycentricRWG_coefficients[:, 6] = 1.
    RWG_to_barycentricRWG_coefficients[:, 7] = 1.
    for i in range(N_RWG):
        # original RWG nodes
        e0 = RWGNumber_edgeVertexes[i, 0]
        e1 = RWGNumber_edgeVertexes[i, 1]
        # now length of RWG
        l_RWG = vertexes_coord_barycentric[e0,:]-vertexes_coord_barycentric[e1,:]
        l = sqrt(sum(l_RWG*l_RWG))
        # we first deal with the first half of the original RWG
        t = RWGNumber_signedTriangles[i, 0]
        nodes_t = divided_triangles_vertexes[t,:].tolist()
        index_start_node = nodes_t.index(e0)
        index_end_node = nodes_t.index(e1)
        index_mid_node = (index_start_node + index_end_node)/2
        if ((index_start_node==0) and (index_end_node==4)) or ((index_start_node==4) and (index_end_node==0)):
            index_mid_node = 5
        # remember:
        # barycentric_RWGNumber_signedTriangles[6*t + 0, :] = [6*t+0, 6*t+1]
        # barycentric_RWGNumber_signedTriangles[6*t + 1, :] = [6*t+1, 6*t+2]
        # barycentric_RWGNumber_signedTriangles[6*t + 2, :] = [6*t+2, 6*t+3]
        # barycentric_RWGNumber_signedTriangles[6*t + 3, :] = [6*t+3, 6*t+4]
        # barycentric_RWGNumber_signedTriangles[6*t + 4, :] = [6*t+4, 6*t+5]
        # barycentric_RWGNumber_signedTriangles[6*t + 5, :] = [6*t+5, 6*t+0]
        # now lengths of inner RWGs
        l_barycentricRWGs = zeros(6, 'd')
        for j in range(6):
            vector = vertexes_coord_barycentric[nodes_t[(j+1)%6]] - vertexes_coord_barycentric[nodes_t[6]]
            l_barycentricRWGs[j] = sqrt(sum(vector*vector))
        coefficients = array([l/6., l/3., l/6., l/3., l/6., l/3.])
        coefficients = coefficients/l_barycentricRWGs
        if index_mid_node == 1:
            # then inner RWGs 1 and 4 (indexes 0 and 3) are aligned with mother RWG and remain zero
            # we only care about inner RWGs 2, 3, 5, 6 (indexes 1, 2, 4, 5)
            # in this case, RWGs 2 and 3 are against mother RWG, and 5 and 6 are in same direction
            RWG_to_barycentricRWG_coefficients[i, :6] = array([0., -1., -1., 0., 1., 1.]) * coefficients
        if index_mid_node == 3:
            # then inner RWGs 3 and 6 (indexes 2 and 5) are aligned with mother RWG and remain zero
            # we only care about inner RWGs 1, 2, 4, 5 (indexes 0, 1, 3, 4)
            # in this case, RWGs 4 and 5 are against mother RWG, and 1 and 2 are in same direction
            RWG_to_barycentricRWG_coefficients[i, :6] = array([1., 1., 0., -1., -1., 0.]) * coefficients
        if index_mid_node == 5:
            # then inner RWGs 2 and 5 (indexes 1 and 4) are aligned with mother RWG and remain zero
            # we only care about inner RWGs 1, 3, 4, 6 (indexes 0, 2, 3, 5)
            # in this case, RWGs 1 and 6 are against mother RWG, and 3 and 4 are in same direction
            RWG_to_barycentricRWG_coefficients[i, :6] = array([-1., 0., 1., 1., 0., -1.]) * coefficients

        # we now deal with the second half of the original RWG
        t = RWGNumber_signedTriangles[i, 1]
        nodes_t = divided_triangles_vertexes[t,:].tolist()
        index_start_node = nodes_t.index(e0)
        index_end_node = nodes_t.index(e1)
        index_mid_node = (index_start_node + index_end_node)/2
        if ((index_start_node==0) and (index_end_node==4)) or ((index_start_node==4) and (index_end_node==0)):
            index_mid_node = 5
        l_barycentricRWGs = zeros(6, 'd')
        for j in range(6):
            vector = vertexes_coord_barycentric[nodes_t[(j+1)%6]] - vertexes_coord_barycentric[nodes_t[6]]
            l_barycentricRWGs[j] = sqrt(sum(vector*vector))
        coefficients = array([l/6., l/3., l/6., l/3., l/6., l/3.])
        coefficients = coefficients/l_barycentricRWGs
        if index_mid_node == 1:
            # then inner RWGs 1 and 4 (indexes 0 and 3) are aligned with mother RWG and remain zero
            # we only care about inner RWGs 2, 3, 5, 6 (indexes 1, 2, 4, 5)
            # in this case, RWGs 2 and 3 are against mother RWG, and 5 and 6 are in same direction
            RWG_to_barycentricRWG_coefficients[i, 8:] = array([0., -1., -1., 0., 1., 1.]) * coefficients * (-1.0)
        if index_mid_node == 3:
            # then inner RWGs 3 and 6 (indexes 2 and 5) are aligned with mother RWG and remain zero
            # we only care about inner RWGs 1, 2, 4, 5 (indexes 0, 1, 3, 4)
            # in this case, RWGs 4 and 5 are against mother RWG, and 1 and 2 are in same direction
            RWG_to_barycentricRWG_coefficients[i, 8:] = array([1., 1., 0., -1., -1., 0.]) * coefficients * (-1.0)
        if index_mid_node == 5:
            # then inner RWGs 2 and 5 (indexes 1 and 4) are aligned with mother RWG and remain zero
            # we only care about inner RWGs 1, 3, 4, 6 (indexes 0, 2, 3, 5)
            # in this case, RWGs 1 and 6 are against mother RWG, and 3 and 4 are in same direction
            RWG_to_barycentricRWG_coefficients[i, 8:] = array([-1., 0., 1., 1., 0., -1.]) * coefficients * (-1.0)

    return RWG_to_barycentricRWG, RWG_to_barycentricRWG_coefficients



if __name__=="__main__":
    path = './geo'
    targetName = 'strip'
    f = 1.12e9
    write_geo(path, targetName, 'lc', c/f/10.0)
    write_geo(path, targetName, 'lx', .05)
    write_geo(path, targetName, 'ly', .05)
    write_geo(path, targetName, 'lz', .05)
    write_geo(path, targetName, 'w', 0.02)
    executeGmsh(path, targetName, 0)
    targetDimensions_scaling_factor = 1.0
    z_offset = 0.0
    vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_1(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)
    #vertexes_coord, triangle_vertexes, triangles_physicalSurface = read_mesh_GMSH_2(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)
    edgeNumber_vertexes, edgeNumber_triangles, triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = edges_computation(triangle_vertexes, vertexes_coord)
    T = len(triangle_adjacentTriangles)
    
    #print "attribution of a surface number to each triangle...",
    #triangles_surfaces = triangles_surfaces_computation(triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction)
    #S = max(triangles_surfaces)+1

    t0 = time.clock()
    triangles_surfaces = reorder_triangle_vertexes(triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction, triangle_vertexes, vertexes_coord)

    print("    checking open and closed surfaces...")
    is_closed_surface, connected_surfaces, potential_closed_surfaces = is_surface_closed(triangles_surfaces, edgeNumber_triangles)
    print("is_closed_surface = " + str(is_closed_surface * 1))
    print("    connected surfaces = " + str(connected_surfaces))
    print("    potential closed surfaces = " + str(potential_closed_surfaces))

    time_reordering_normals = time.clock()-t0
    print("time = " + str(time_reordering_normals) + " seconds")
    #triangles_centroids = triangles_centroids_computation(vertexes_coord, triangle_vertexes)
    #triangles_areas, triangles_normals = triangles_areas_normals_computation(vertexes_coord, triangle_vertexes, triangles_surfaces)
    #write_normals(os.path.join(path, "normals.pos"), triangles_centroids, triangles_normals, triangles_surfaces, -1)
    
    RWGNumber_signedTriangles, RWGNumber_edgeVertexes, N_edges, N_RWG = RWGNumber_signedTriangles_computation(edgeNumber_triangles, edgeNumber_vertexes, triangles_surfaces, is_closed_surface, triangle_vertexes, vertexes_coord)
    #print RWGNumber_signedTriangles
    RWGNumber_oppVertexes = RWGNumber_oppVertexes_computation(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, triangle_vertexes)
    print("    Number of edges = " + str(N_edges))
    print("    Number of RWG = " + str(N_RWG))


    divided_triangles_vertexes, MAX_V = divide_triangles(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, triangle_vertexes, vertexes_coord)
    vertexes_coord_barycentric, barycentric_triangles_vertexes = create_barycentric_triangles(divided_triangles_vertexes, vertexes_coord, MAX_V)
    print "number of original triangles = ", triangle_vertexes.shape[0]
    print "number of barycentric triangles = ", barycentric_triangles_vertexes.shape[0]
    print vertexes_coord_barycentric.shape

    barycentric_RWGNumber_signedTriangles, barycentric_RWGNumber_edgeVertexes, barycentric_RWGNumber_oppVertexes = create_barycentric_RWGs(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, divided_triangles_vertexes, barycentric_triangles_vertexes)
    T = divided_triangles_vertexes.shape[0]
    N_RWG_bary = barycentric_RWGNumber_signedTriangles.shape[0]
    print "N_RWG_bary =", N_RWG_bary
    RWG_to_barycentricRWG, RWG_to_barycentricRWG_coefficients = create_RWG_to_barycentricRWG(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, divided_triangles_vertexes, vertexes_coord_barycentric)
    print RWG_to_barycentricRWG_coefficients
