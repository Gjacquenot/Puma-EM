import os, sys
from math import pi
from scipy import zeros, array, arange, dot
from scipy import sparse, linalg, cos, sin, conj, log10, real, sum, imag
from scipy.sparse.linalg import gmres, bicgstab, lgmres
from meshClass import MeshClass
from PyGmsh import executeGmsh, write_geo
from Z_MoM import Z_CFIE_MoM, Z_EH_J_MoM
from V_EH import computeV_EH, G_EJ_G_HJ, V_EH_dipole_alternative, V_EH_plane
from EM_constants import *
from MoMPostProcessing import *
try:
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    pass
from mesh_functions_seb import *


def computeCurrentsVisualization(w, target_mesh, ZI):
    if (N_RWG<1e4):
        nbTimeSteps = 48
    else:
        nbTimeSteps = 1
    J_centroids_triangles = JMCentroidsTriangles(ZI, target_mesh)
    norm_J_centroids_triangles = normJMCentroidsTriangles(J_centroids_triangles, w, nbTimeSteps)
    write_VectorFieldTrianglesCentroids(os.path.join("./geo", target_mesh.targetName) + '.J_centroids_triangles.pos', real(J_centroids_triangles), target_mesh)
    write_ScalarFieldTrianglesCentroids(os.path.join("./geo", target_mesh.targetName) + '.norm_J_centroids_triangles.pos', norm_J_centroids_triangles, target_mesh)

class Target_MoM:
    def __init__(self, CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION):
        print "Target_MoM instanciation..."
        signSurfObs, signSurfSrc = 1.0, 1.0 # no dielectric target here
        target_mesh.RWGNumber_M_CURRENT_OK *= 0
        self.Z_CFIE_J, self.Z_CFIE_M = Z_CFIE_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
        self.iter_counter = 0
    # functions
    def matvec(self, x):
        self.iter_counter += 1
        return dot(self.Z_CFIE_J, x)

    def V_EH_computation(self, CFIE, target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION):
        self.V_EH = computeV_EH(target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION, 'F')
        self.V_CFIE = self.V_EH[:,0] * CFIE[0]
        for i in range(1, 4):
            self.V_CFIE += self.V_EH[:,i] * CFIE[i] * target_mesh.RWGNumber_CFIE_OK

    def compute_Y_CFIE(self):
        self.Y_CFIE = linalg.inv(self.Z_CFIE_J)

    def solveByInversion(self):
       print "Matrix inversion..."
       self.compute_Y_CFIE()
       self.I_CFIE = dot(self.Y_CFIE, self.V_CFIE)

    def solveByLUdecomposition(self):
        print "LU decomposition and solution..."
        t0 = time.clock()
        lu, piv = linalg.lu_factor(self.Z_CFIE_J)
        self.I_CFIE = linalg.lu_solve((lu, piv), self.V_CFIE)
        print "Done. time =", time.clock() - t0, "seconds"


def itercount(residual):
    global count
    count = count + 1

if __name__=="__main__":
    path = './geo'
    targetName = 'cube2'
    # first resonances for the sphere: f = c * Z/(2*pi*a), where Z is a zero of J Bessel function
    # Z = 4.493409375, 5.763459195, 6.987932, 8.18256145, 9.35581211, 10.5128354
    # respectively for orders 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5
    # Frob-EFIE convergence is difficult or impossible at the corresponding frequencies,
    # especially for order 5.5, a = 0.3, for which f = 1487993627.3926289, tol = 1e-3
    # However, Frob-CFIE convergence is more than OK: it is guaranteed
    f = .9e9
    fileName = targetName
    write_geo(path, fileName, 'lc', c/f/10.)
    write_geo(path, fileName, 'lx', 0.075)
    write_geo(path, fileName, 'ly', 0.07)
    write_geo(path, fileName, 'lz', 0.07)
    executeGmsh(path, targetName, 0)
    z_offset = 0.0
    targetDimensions_scaling_factor = 1.0
    languageForMeshConstruction = "C++"
    meshFormat = 'GMSH' 
    meshFileTermination = '.msh'
    target_mesh = MeshClass(path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction, meshFormat, meshFileTermination)
    target_mesh.constructFromGmshFile()
    N_RWG = target_mesh.N_RWG
    print("average RWG length = " + str(target_mesh.average_RWG_length) + "m = lambda /" + str((c/f)/target_mesh.average_RWG_length))
    print 
    w = 2. * pi * f
    eps_r = 1.
    mu_r = 1.
    TDS_APPROX = 0
    Z_s = 0.0
    J_dip = array([1.0, 1.0, 0.], 'D')
    r_dip = array([0.1, 0.1, 2.0], 'd')
    list_of_test_edges_numbers = arange(N_RWG,dtype='i')
    list_of_src_edges_numbers = arange(N_RWG,dtype='i')
    MOM_FULL_PRECISION = 1
    EXCITATION = 'dipole'
    CHOICE = "CFIE testing"
    #CHOICE = "fields verification"
    #CHOICE = "dielectric target"
    V_RWG = zeros(N_RWG, 'D')

    if CHOICE=="CFIE testing":
        #for coeff in [1.0, .8, 0.5, 0.2, 0.0]:
        for coeff in [1.0]:
            #CFIE = array([coeff, coeff, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            #CFIE = array([coeff, 0, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            CFIE = array([coeff* 1.0/sqrt(mu_0/eps_0), 0, -(1.0 - coeff), -(1.0 - coeff)], 'D')
            print "CFIE =", CFIE
            target_MoM = Target_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
            # excitation computation
            target_MoM.V_EH_computation(CFIE, target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION)
            V_RWG = target_MoM.V_CFIE
            #target_MoM.solveByInversion()
            target_MoM.solveByLUdecomposition()
            print "inverted MoM RCS =", sum(target_MoM.I_CFIE*target_MoM.V_EH[:,0])
            #computeCurrentsVisualization(w, target_mesh, target_MoM.I_CFIE)
            # now we try the iterative method
            count = 0
            I_CFIE_bicgstab = bicgstab(target_MoM.Z_CFIE_J, target_MoM.V_CFIE, x0=None, tol=1.0e-05, maxiter=1000, xtype=None, callback=itercount)
            print "bicgstab MoM RCS =", sum(I_CFIE_bicgstab[0]*target_MoM.V_EH[:,0]), "# of iterations =", count

    # now construction of the barycentric mesh: the classical way
    target_mesh_bary = MeshClass(path, targetName, targetDimensions_scaling_factor, z_offset, languageForMeshConstruction, meshFormat, meshFileTermination)
    divided_triangles_vertexes, MAX_V = divide_triangles(target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.triangle_vertexes, target_mesh.vertexes_coord)
    target_mesh_bary.vertexes_coord, target_mesh_bary.triangle_vertexes = create_barycentric_triangles(divided_triangles_vertexes, target_mesh.vertexes_coord, MAX_V)

    #t0 = time.clock()
    edgeNumber_vertexes, edgeNumber_triangles, triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = edges_computation(target_mesh_bary.triangle_vertexes, target_mesh_bary.vertexes_coord)
    #print is_triangle_adjacentTriangles_via_junction
    #time_edges_classification = time.clock()-t0
    #print("  edges classification cumulated time = " + str(time_edges_classification) + " seconds")

    #print("  reordering triangles for normals coherency...")
    #sys.stdout.flush()
    #t0 = time.clock()
    target_mesh_bary.triangles_surfaces = reorder_triangle_vertexes(triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction, target_mesh_bary.triangle_vertexes, target_mesh_bary.vertexes_coord)
    target_mesh_bary.S = max(target_mesh_bary.triangles_surfaces)+1
    #time_reordering_normals = time.clock()-t0
    #print("  cumulated time = " + str(time_reordering_normals) + " seconds")

    #print("  checking the closed and open surfaces...")
    #sys.stdout.flush()
    #t0 = time.clock()
    target_mesh_bary.IS_CLOSED_SURFACE, target_mesh_bary.connected_surfaces, target_mesh_bary.potential_closed_surfaces = is_surface_closed(target_mesh_bary.triangles_surfaces, edgeNumber_triangles)
    print("  test of the closed surfaces : " + str(target_mesh_bary.IS_CLOSED_SURFACE))
    print("  connected surfaces : " + str(target_mesh_bary.connected_surfaces))
    print("  potential closed surfaces : " + str(target_mesh_bary.potential_closed_surfaces))

    #print("  computing the effective RWG functions and their opposite vertexes...")
    #sys.stdout.flush()
    #t0 = time.clock()
    target_mesh_bary.RWGNumber_signedTriangles, target_mesh_bary.RWGNumber_edgeVertexes, target_mesh_bary.N_edges, target_mesh_bary.N_RWG = RWGNumber_signedTriangles_computation(edgeNumber_triangles, edgeNumber_vertexes, target_mesh_bary.triangles_surfaces, target_mesh_bary.IS_CLOSED_SURFACE, target_mesh_bary.triangle_vertexes, target_mesh_bary.vertexes_coord)
    del edgeNumber_vertexes
    target_mesh_bary.RWGNumber_oppVertexes = RWGNumber_oppVertexes_computation(target_mesh_bary.RWGNumber_signedTriangles, target_mesh_bary.RWGNumber_edgeVertexes, target_mesh_bary.triangle_vertexes)
    ## memory-economic way for computing average_RWG_length
    #time_effective_RWG_functions_computation =  time.clock() - t0
    #print("  effective RWG functions computation cumulated time = " + str(time_effective_RWG_functions_computation))
    print("  Number of edges = " + str(target_mesh_bary.N_edges))
    print("  Number of RWG = " + str(target_mesh_bary.N_RWG))
    #sys.stdout.flush()
    target_mesh_bary.compute_RWG_CFIE_OK()
    if target_mesh_bary.N_RWG<1e4:
        stride = 1
    else:
        stride = target_mesh_bary.N_RWG/100
    target_mesh_bary.average_RWG_length = compute_RWG_meanEdgeLength(target_mesh_bary.vertexes_coord, target_mesh_bary.RWGNumber_edgeVertexes, stride)
    print("average barycentric RWG length = " + str(target_mesh_bary.average_RWG_length) + "m = lambda /" + str((c/f)/target_mesh_bary.average_RWG_length))


    list_of_test_edges_numbers = arange(target_mesh_bary.N_RWG,dtype='i')
    list_of_src_edges_numbers = arange(target_mesh_bary.N_RWG,dtype='i')
    MOM_FULL_PRECISION = 1
    EXCITATION = 'dipole'
    CHOICE = "CFIE testing"
    #CHOICE = "fields verification"
    #CHOICE = "dielectric target"
    if CHOICE=="CFIE testing":
        #for coeff in [1.0, .8, 0.5, 0.2, 0.0]:
        for coeff in [1.0]:
            #CFIE = array([coeff, coeff, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            #CFIE = array([coeff, 0, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            CFIE = array([coeff* 1.0/sqrt(mu_0/eps_0), 0, -(1.0 - coeff), -(1.0 - coeff)], 'D')
            print "CFIE =", CFIE
            target_MoM = Target_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh_bary, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
            # excitation computation
            target_MoM.V_EH_computation(CFIE, target_mesh_bary, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION)
            #target_MoM.solveByInversion()
            target_MoM.solveByLUdecomposition()
            print "inverted MoM RCS =", sum(target_MoM.I_CFIE*target_MoM.V_EH[:,0])
            #computeCurrentsVisualization(w, target_mesh, target_MoM.I_CFIE)
            # now we try the iterative method
            count = 0
            I_CFIE_bicgstab = bicgstab(target_MoM.Z_CFIE_J, target_MoM.V_CFIE, x0=None, tol=1.0e-05, maxiter=1000, xtype=None, callback=itercount)
            print "bicgstab MoM RCS =", sum(I_CFIE_bicgstab[0]*target_MoM.V_EH[:,0]), "# of iterations =", count


    # now construction of the barycentric mesh: the new way based on the original mesh
    target_mesh_bary.RWGNumber_signedTriangles, target_mesh_bary.RWGNumber_edgeVertexes, target_mesh_bary.RWGNumber_oppVertexes = create_barycentric_RWGs(target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, divided_triangles_vertexes, target_mesh_bary.triangle_vertexes)
    target_mesh_bary.N_RWG = target_mesh_bary.RWGNumber_signedTriangles.shape[0]
    #print("  Number of edges = " + str(target_mesh_bary.N_edges))
    print("  Number of RWG = " + str(target_mesh_bary.N_RWG))
    sys.stdout.flush()
    if target_mesh_bary.N_RWG<1e4:
        stride = 1
    else:
        stride = target_mesh_bary.N_RWG/100
    target_mesh_bary.average_RWG_length = compute_RWG_meanEdgeLength(target_mesh_bary.vertexes_coord, target_mesh_bary.RWGNumber_edgeVertexes, stride)
    print("average barycentric RWG length = " + str(target_mesh_bary.average_RWG_length) + "m = lambda /" + str((c/f)/target_mesh_bary.average_RWG_length))
    T = divided_triangles_vertexes.shape[0]
    target_mesh_bary.RWGNumber_CFIE_OK, target_mesh_bary.RWGNumber_M_CURRENT_OK = compute_barycentric_RWG_CFIE_OK(target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, 6*T)
    print "T_bary =", 6*T
    list_of_test_edges_numbers = arange(target_mesh_bary.N_RWG,dtype='i')
    list_of_src_edges_numbers = arange(target_mesh_bary.N_RWG,dtype='i')
    MOM_FULL_PRECISION = 1
    EXCITATION = 'dipole'
    CHOICE = "CFIE testing"
    #CHOICE = "fields verification"
    #CHOICE = "dielectric target"
    V_RWG_bary = zeros(target_mesh_bary.N_RWG, 'D')
    if CHOICE=="CFIE testing":
        #for coeff in [1.0, .8, 0.5, 0.2, 0.0]:
        for coeff in [1.0]:
            #CFIE = array([coeff, coeff, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            #CFIE = array([coeff, 0, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            CFIE = array([coeff* 1.0/sqrt(mu_0/eps_0), 0, -(1.0 - coeff), -(1.0 - coeff)], 'D')
            print "CFIE =", CFIE
            target_MoM = Target_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh_bary, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
            # excitation computation
            target_MoM.V_EH_computation(CFIE, target_mesh_bary, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION)
            V_RWG_bary = target_MoM.V_CFIE
            #target_MoM.solveByInversion()
            target_MoM.solveByLUdecomposition()
            print "inverted MoM RCS =", sum(target_MoM.I_CFIE*target_MoM.V_EH[:,0])
            #computeCurrentsVisualization(w, target_mesh, target_MoM.I_CFIE)
            # now we try the iterative method
            count = 0
            I_CFIE_bicgstab = bicgstab(target_MoM.Z_CFIE_J, target_MoM.V_CFIE, x0=None, tol=1.0e-05, maxiter=1000, xtype=None, callback=itercount)
            print "bicgstab MoM RCS =", sum(I_CFIE_bicgstab[0]*target_MoM.V_EH[:,0]), "# of iterations =", count


    RWG_to_barycentricRWG, RWG_to_barycentricRWG_coefficients = create_RWG_to_barycentricRWG(target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, divided_triangles_vertexes, target_mesh_bary.vertexes_coord)

    N_RWG = V_RWG.shape[0]
    N_RWG_bary = V_RWG_bary.shape[0]
    V_RWG_2 = zeros(N_RWG, 'D')
    for i in range(N_RWG):
        for j in range(14):
            index = RWG_to_barycentricRWG[i, j]
            V_RWG_2[i] += RWG_to_barycentricRWG_coefficients[i, j]* V_RWG_bary[index]

    #for i in range(N_RWG):
        #print V_RWG[i], V_RWG_2[i], real(V_RWG_2[i]-V_RWG[i])/real(V_RWG[i]), imag(V_RWG_2[i]-V_RWG[i])/imag(V_RWG[i])


