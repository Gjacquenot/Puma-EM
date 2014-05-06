from math import pi
from scipy import real, imag, exp
from mesh_functions_seb import *

def JMCentroidsTriangles(ICoeffs, target_mesh):
    """given a mesh and MoM current coefficients, this function
    aims at calculating the current at each centroid"""
    J_M_centroids = zeros((target_mesh.T, 3), 'D')
    triangles_centroids = triangles_centroids_computation(target_mesh.vertexes_coord, target_mesh.triangle_vertexes)
    triangles_areas, triangles_normals = triangles_areas_normals_computation(target_mesh.vertexes_coord, target_mesh.triangle_vertexes, target_mesh.triangles_surfaces)
    RWGNumber_edgeLength = compute_RWGNumber_edgeLength(target_mesh.vertexes_coord, target_mesh.RWGNumber_edgeVertexes)
    RWG_number = 0
    for elem in target_mesh.RWGNumber_signedTriangles:
        for i in range(len(elem)):
            t = abs(elem[i])
            rCentroid = triangles_centroids[t]
            Area = triangles_areas[t]
            signRWG = (-1)**(i)
            rp = target_mesh.vertexes_coord[target_mesh.RWGNumber_oppVertexes[RWG_number, i]]
            lp = RWGNumber_edgeLength[RWG_number]
            f_triangle_i_edge_p = lp/(2.0*Area)*signRWG*(rCentroid-rp)
            J_M_centroids[t] += ICoeffs[RWG_number] * f_triangle_i_edge_p
        RWG_number += 1
    return J_M_centroids
    

def normJMCentroidsTriangles(J_M_centroids, w, nbTimeSteps):
    timeValues = arange(0, 2*pi/w, 2*pi/w/nbTimeSteps)
    norm_J_M_centroids = zeros((J_M_centroids.shape[0], nbTimeSteps), 'f')
    # the norm is computed at time t -> J(t) = real(J * exp(jwt))
    # if t=0, J(t) = real(J)
    for k in range(nbTimeSteps):
        JTimeDomain = real(J_M_centroids * exp(1.j * w * timeValues[k]))
        norm_J_M_centroids[:, k] = sqrt( sum(JTimeDomain**2, 1) )
    return norm_J_M_centroids

def write_VectorFieldTrianglesCentroids(name, VectorField, target_mesh):
    """function that writes a given vector field on the triangles centroids."""
    triangles_centroids = triangles_centroids_computation(target_mesh.vertexes_coord, target_mesh.triangle_vertexes)
    f = open(name, 'w')
    for k in range(target_mesh.T):
        string_to_write = ''
        for i in range(3):
            string_to_write += str(real(VectorField[k][i])) + ' ' + str(imag(VectorField[k][i])) + ' '
        for i in range(3):
            string_to_write += str(triangles_centroids[k][i]) + ' '
        string_to_write += '\n'
        f.write(string_to_write)
    f.close()

def write_VectorFieldTrianglesCentroidsGMSH(name, VectorField, target_mesh):
    """function that writes a given vector field on a given surface
    to a file readable by GMSH for viewing."""
    triangles_centroids = triangles_centroids_computation(target_mesh.vertexes_coord, target_mesh.triangle_vertexes)
    f = open(name, 'w')
    f.write('View "surface currents of surfaces xxx" {\n')
    for k in range(target_mesh.T):
        write_condition = 1
        if write_condition:
            string_to_write = 'VP(' + str(triangles_centroids[k].tolist())[1:-1] + ')'
            string_to_write += '{' + str(VectorField[k].tolist())[1:-1] + '};\n'
            f.write(string_to_write)
    f.write('};\n')
    f.close()

def write_ScalarFieldTrianglesCentroidsGMSH(name, ScalarField, target_mesh):
    """function that writes a given scalar field on a given surface
    to a file readable by GMSH for viewing."""
    f = open(name, 'w')
    f.write('View "surface currents of surfaces xxx" {\n')
    Nc = ScalarField.shape[1]
    for i in range(target_mesh.T):
        string_to_write = 'ST('
        r0 = target_mesh.vertexes_coord[target_mesh.triangle_vertexes[i, 0]]
        r1 = target_mesh.vertexes_coord[target_mesh.triangle_vertexes[i, 1]]
        r2 = target_mesh.vertexes_coord[target_mesh.triangle_vertexes[i, 2]]
        string_to_write += str(r0.tolist())[1:-1]
        string_to_write += ', ' + str(r1.tolist())[1:-1]
        string_to_write += ', ' + str(r2.tolist())[1:-1] + ')'
        string_to_write += '{'
        for j in range(Nc-1):
            string_to_write += str(ScalarField[i, j]) + ', ' + str(ScalarField[i, j]) + ', ' +  str(ScalarField[i, j]) + ', '
        string_to_write += str(ScalarField[i, Nc-1]) + ', ' + str(ScalarField[i, Nc-1]) + ', ' +  str(ScalarField[i, Nc-1])
        string_to_write += '};\n'
        f.write(string_to_write)
    f.write('};\n')
    f.close()

