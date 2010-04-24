import os
from math import pi
from scipy import zeros, array, arange, dot
from scipy import linalg, cos, sin, conj, log10, real, sum, imag
from scipy.linalg import gmres, bicgstab
from meshClass import MeshClass
from PyGmsh import executeGmsh, write_geo
from Z_MoM import Z_MoM
from V_EH import computeV_EH, G_EJ_G_HJ, V_EH_dipole_alternative
from EM_constants import *
from MoMPostProcessing import *

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
        self.Z_CFIE_J, self.Z_CFIE_M = Z_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
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

class dielectricTarget_MoM:
    def __init__(self, CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r_out, mu_r_out, eps_r_in, mu_r_in, MOM_FULL_PRECISION):
        self.numberOfMedia = sum(target_mesh.IS_CLOSED_SURFACE) + 1
        print "MOM.py: number of possible media =", self.numberOfMedia
        print "Target_MoM instanciation..."
        TDS_APPROX, Z_s = 0, 0.0 + 0.0j
        N_J, N_M = len(list_of_test_edges_numbers), len(list_of_test_edges_numbers)
        self.Z_CFIE = zeros((N_J + N_M, N_J + N_M), 'D')
        print "dielectricTarget_MoM: computing the outside interactions..."
        signSurfObs, signSurfSrc = 1.0, 1.0
        self.Z_CFIE[:N_J, :N_J], self.Z_CFIE[:N_J,N_J:N_J + N_M] = Z_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r_out, mu_r_out, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
        print "dielectricTarget_MoM: computing the inside interactions..."
        signSurfObs, signSurfSrc = -1.0, -1.0
        self.Z_CFIE[N_J:N_J + N_M,:N_J], self.Z_CFIE[N_J:N_J + N_M, N_J:N_J + N_M] = Z_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r_in, mu_r_in, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
        self.iter_counter = 0
    # functions
    def matvec(self, x):
        self.iter_counter += 1
        return dot(self.Z_CFIE, x)

    def V_EH_computation(self, CFIE, target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION):
        V_EH_tmp = computeV_EH(target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION, 'F')
        self.V_EH = zeros((len(list_of_test_edges_numbers) + len(list_of_test_edges_numbers), 4), 'D')
        self.V_EH[:len(list_of_test_edges_numbers)] = V_EH_tmp
        self.V_CFIE = self.V_EH[:,0] * CFIE[0]
        for i in range(1, 4):
            self.V_CFIE += self.V_EH[:,i] * CFIE[i]
    
    def compute_Y_CFIE(self):
        self.Y_CFIE = linalg.inv(self.Z_CFIE)

    def solveByInversion(self):
        print "Matrix inversion..."
        self.compute_Y_CFIE()
        self.I_CFIE = dot(self.Y_CFIE, self.V_CFIE)

if __name__=="__main__":
    path = './geo'
    targetName = 'sphere'
    # first resonances for the sphere: f = c * Z/(2*pi*a), where Z is a zero of J Bessel function
    # Z = 4.493409375, 5.763459195, 6.987932, 8.18256145, 9.35581211, 10.5128354
    # respectively for orders 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5
    # Frob-EFIE convergence is difficult or impossible at the corresponding frequencies,
    # especially for order 5.5, a = 0.3, for which f = 1487993627.3926289, tol = 1e-3
    # However, Frob-CFIE convergence is more than OK: it is guaranteed
    f = 2.12e9
    fileName = targetName
    write_geo(path, fileName, 'lc', c/f/15.)
    write_geo(path, fileName, 'lx', 0.07)
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

    w = 2. * pi * f
    eps_r = 1.
    mu_r = 1.
    TDS_APPROX = 0
    Z_s = 0.0
    J_dip = array([.0, 1.0, 0.], 'D')
    r_dip = array([0.0, 0.0, 0.2714], 'd')
    list_of_test_edges_numbers = arange(N_RWG,dtype='i')
    list_of_src_edges_numbers = arange(N_RWG,dtype='i')
    MOM_FULL_PRECISION = 1
    EXCITATION = 'dipole'
    CHOICE = "CFIE testing"
    #CHOICE = "fields verification"
    #CHOICE = "dielectric target"
    if CHOICE=="CFIE testing":
        for coeff in [1.0, 0.8, 0.5, 0.2, 0.0]:
        #for coeff in [0.2]:
            CFIE = array([coeff, 0.0, 0.0, -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            print
            print "CFIE =", CFIE
            target_MoM = Target_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
            # excitation computation
            target_MoM.V_EH_computation(CFIE, target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION)
            target_MoM.solveByInversion()
            #computeCurrentsVisualization(w, target_mesh, target_MoM.I_CFIE)
            # now we try the iterative method
            I_CFIE_bicgstab = bicgstab(target_MoM, target_MoM.V_CFIE, x0=None, tol=1.0e-05, maxiter=1000, xtype=None)
            print "bicgstab MoM RCS =", sum(I_CFIE_bicgstab[0]*target_MoM.V_EH[:,0]), "# of iterations =", target_MoM.iter_counter

    if CHOICE=="fields verification":
        coeff = 0.2
        CFIE = array([coeff, 0.0, 0.0, -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
        target_MoM = Target_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
        # excitation computation
        target_MoM.V_EH_computation(CFIE, target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION)
        target_MoM.solveByInversion()

        positions = arange(0,0.5,0.01).astype('d')
        H_x_obs = zeros(len(positions), 'D')
        index = 0
        for pos in positions:
            r_obs = array([0., pos, 0.], 'd')
            J_obs = array([1., 0., 0.], 'D')
            V_EH_obs = V_EH_dipole_alternative(J_obs, r_obs, list_of_test_edges_numbers, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r)
            H_x_obs[index] = dot(target_MoM.I_CFIE, V_EH_obs[:,2])
            index += 1
        
        from pylab import plot, rc, subplot, xlabel, ylabel, legend, xticks, yticks, grid, gca, setp, show
        #rc('text', usetex=True)
        FontSize=18
        LineWidth=1
        plot(positions, real(H_x_obs),'b',positions,imag(H_x_obs),'r--')
        xlabel(r'distance from center',fontsize=FontSize+2)
        ylabel(r'$E_x$ field',fontsize=FontSize+2)
        legend([r'real$E_x$', r'imag$E_x$'])
        xticks(fontsize=FontSize)
        yticks(fontsize=FontSize)
        grid(True)
        show()

    if CHOICE=="dielectric target":
        coeff = 0.2
        eps_r_out, mu_r_out = eps_r, mu_r
        eps_r_in, mu_r_in = 3.0 - 2.0j, 1.0 - 0.0j
        #CFIE = array([coeff, 0.0, 0., -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
        TE, NE, TH, NH = 1., 0., 1., 1.
        CFIE = array([TE * coeff, NE * coeff, -TH * (1.0 - coeff) * sqrt(mu_0/eps_0), -NH * (1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
        target_MoM = dielectricTarget_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r_out, mu_r_out, eps_r_in, mu_r_in, MOM_FULL_PRECISION)
        N_J, N_M = len(list_of_test_edges_numbers), len(list_of_test_edges_numbers)
        
        
        J_dip = array([0.0, 1.0, 0.], 'D') * sqrt((c/f)/sqrt(mu_0/eps_0) * 3.0/pi)
        r_dip = array([0.0, 0.0, 0.2714], 'd')

        target_MoM.V_EH_computation(CFIE, target_mesh, J_dip, r_dip, w, eps_r_out, mu_r_out, list_of_test_edges_numbers, EXCITATION)
        target_MoM.solveByInversion()
        print sum(target_MoM.I_CFIE)

        ##I_CFIE_iter = gmres(target_MoM, target_MoM.V_CFIE, x0=None, tol=1.0e-03, maxiter=1000, xtype=None, restrt=100)[0]
        
        ## computation of the scattered fields
        J_obs = array([1.0, 0.0, 0.0], 'D')
        r_obs = array([0.0, .12, 0.0], 'd')
        V_EH = computeV_EH(target_mesh, J_obs, r_obs, w, eps_r_out, mu_r_out, list_of_test_edges_numbers, 'dipole', 'F')
        
        # incoming field computation
        G_EJ, G_HJ = G_EJ_G_HJ(r_dip, r_obs, eps_r_out, mu_r_out, w/c)
        E_inc, H_inc = dot(G_EJ, J_obs), dot(G_HJ, J_obs)
        print "incoming E_field at observation point =", E_inc
        
        print "MoM RCS =", sum(-target_MoM.I_CFIE[:N_J]*V_EH[:,0] + target_MoM.I_CFIE[N_J:N_J + N_M]*V_EH[:,2])
        print "total =", sum(-target_MoM.I_CFIE[:N_J]*V_EH[:,0] + target_MoM.I_CFIE[N_J:N_J + N_M]*V_EH[:,2]) + E_inc[0]
        #print "bicgstab dielectric MoM RCS =", sum(-I_CFIE_iter[:N_J]*V_EH[:,0] + I_CFIE_iter[N_J:N_J + N_M]*V_EH[:,2]) + E_inc[0], "# of iterations =", target_MoM.iter_counter



