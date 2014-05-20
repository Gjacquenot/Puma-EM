import os
from math import pi
from scipy import zeros, array, arange, dot
from scipy import sparse, linalg, cos, sin, conj, log10, real, sum, imag
from scipy.sparse.linalg import gmres, bicgstab
from meshClass import MeshClass
from PyGmsh import executeGmsh, write_geo
from Z_MoM import Z_MoM
from V_EH import computeV_EH, G_EJ_G_HJ, V_EH_dipole_alternative, V_EH_plane
from EM_constants import *
from MoMPostProcessing import *
try:
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    pass

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

    def solveByLUdecomposition(self):
        print "LU decomposition and solution..."
        t0 = time.clock()
        lu, piv = linalg.lu_factor(self.Z_CFIE_J)
        self.I_CFIE = linalg.lu_solve((lu, piv), self.V_CFIE)
        print "Done. time =", time.clock() - t0, "seconds"

class dielectricTarget_MoM:
    def __init__(self, CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r_out, mu_r_out, eps_r_in, mu_r_in, MOM_FULL_PRECISION, IS_PMCHWT):
        self.numberOfMedia = sum(target_mesh.IS_CLOSED_SURFACE) + 1
        print "MOM.py: number of possible media =", self.numberOfMedia
        print "Target_MoM instanciation..."
        t0 = time.clock()
        TDS_APPROX, Z_s = 0, 0.0 + 0.0j
        N_J, N_M = len(list_of_test_edges_numbers), len(list_of_test_edges_numbers)
        self.Z = zeros((N_J + N_M, N_J + N_M), 'D')
        print "dielectricTarget_MoM: computing the outside interactions..."
        signSurfObs, signSurfSrc = 1.0, 1.0

        if not (IS_PMCHWT == 1): # we have CFIE
            self.Z[:N_J, :N_J], self.Z[:N_J,N_J:N_J + N_M] = Z_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r_out, mu_r_out, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)

        else: # we have PMCHWT
            CFIE_for_PMCHWT = array([1.0, 0., 0., 0.], 'D')
            Z_EJ, Z_EM = Z_MoM(CFIE_for_PMCHWT, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r_out, mu_r_out, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
            self.Z[:N_J, :N_J], self.Z[:N_J,N_J:N_J + N_M] = Z_EJ, Z_EM
            self.Z[N_J:N_J+N_M, :N_J] = -Z_EM # Z_HJ = -Z_EM
            self.Z[N_J:N_J+N_M,N_J:N_J + N_M] = (eps_0*eps_r_out / (mu_0 * mu_r_out)) * Z_EJ # Z_HM = eps/mu * Z_EJ

        print "dielectricTarget_MoM: computing the inside interactions..."
        signSurfObs, signSurfSrc = -1.0, -1.0

        if not (IS_PMCHWT == 1): # we have CFIE
            self.Z[N_J:N_J + N_M,:N_J], self.Z[N_J:N_J + N_M, N_J:N_J + N_M] = Z_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r_in, mu_r_in, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)

        else: # we have PMCHWT
            CFIE_for_PMCHWT = array([1.0, 0., 0., 0.], 'D')
            Z_EJ, Z_EM = Z_MoM(CFIE_for_PMCHWT, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_M_CURRENT_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r_in, mu_r_in, signSurfObs, signSurfSrc, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
            self.Z[:N_J, :N_J] += Z_EJ
            self.Z[:N_J,N_J:N_J + N_M] += Z_EM
            self.Z[N_J:N_J+N_M, :N_J] -= Z_EM # Z_HJ = -Z_EM
            self.Z[N_J:N_J+N_M,N_J:N_J + N_M] += (eps_0*eps_r_in / (mu_0 * mu_r_in)) * Z_EJ # Z_HM = eps/mu * Z_EJ

        print "Done. time =", time.clock() - t0, "seconds"
        self.iter_counter = 0
    # functions
    def matvec(self, x):
        self.iter_counter += 1
        return dot(self.Z, x)

    def V_EH_computation(self, CFIE, target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION, IS_PMCHWT):
        V_EH_tmp = computeV_EH(target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION, 'F')
        N_test = len(list_of_test_edges_numbers)
        self.V_EH = zeros((2*N_test, 4), 'D')
        self.V_EH[:N_test] = V_EH_tmp
        if not (IS_PMCHWT == 1): # we have CFIE
            self.V = self.V_EH[:,0] * CFIE[0]
            for i in range(1, 4):
                self.V += self.V_EH[:,i] * CFIE[i]
        else:
            self.V = self.V_EH[:,0]
            self.V[N_test:2*N_test] = self.V_EH[:N_test,2]

    def V_EH_plane_computation(self, CFIE, target_mesh, E_0, k_hat, r_ref, w, eps_r, mu_r, list_of_test_edges_numbers, IS_PMCHWT):
        V_EH_tmp = V_EH_plane(E_0, k_hat, r_ref, list_of_test_edges_numbers, target_mesh.RWGNumber_CFIE_OK, target_mesh.RWGNumber_signedTriangles, target_mesh.RWGNumber_edgeVertexes, target_mesh.RWGNumber_oppVertexes, target_mesh.vertexes_coord, w, eps_r, mu_r)
        N_test = len(list_of_test_edges_numbers)
        self.V_EH = zeros((2*N_test, 4), 'D')
        self.V_EH[:N_test] = V_EH_tmp
        if not (IS_PMCHWT == 1): # we have CFIE
            self.V = self.V_EH[:,0] * CFIE[0]
            for i in range(1, 4):
                self.V += self.V_EH[:,i] * CFIE[i]
        else:
            self.V = self.V_EH[:,0]
            self.V[N_test:2*N_test] = self.V_EH[:N_test,2]

    def compute_Y(self):
        self.Y = linalg.inv(self.Z)

    def solveByInversion(self):
        print "Matrix inversion..."
        t0 = time.clock()
        self.compute_Y()
        print "Done. time =", time.clock() - t0, "seconds"
        self.I = dot(self.Y, self.V)

    def solveByLUdecomposition(self):
        print "LU decomposition and solution..."
        t0 = time.clock()
        lu, piv = linalg.lu_factor(self.Z)
        self.I = linalg.lu_solve((lu, piv), self.V)
        print "Done. time =", time.clock() - t0, "seconds"

def itercount(residual):
    global count
    count = count + 1

if __name__=="__main__":
    path = './geo'
    targetName = 'sphere2'
    # first resonances for the sphere: f = c * Z/(2*pi*a), where Z is a zero of J Bessel function
    # Z = 4.493409375, 5.763459195, 6.987932, 8.18256145, 9.35581211, 10.5128354
    # respectively for orders 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5
    # Frob-EFIE convergence is difficult or impossible at the corresponding frequencies,
    # especially for order 5.5, a = 0.3, for which f = 1487993627.3926289, tol = 1e-3
    # However, Frob-CFIE convergence is more than OK: it is guaranteed
    f = 1.5e9
    fileName = targetName
    write_geo(path, fileName, 'lc', c/f/14.)
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

    w = 2. * pi * f
    eps_r = 1.
    mu_r = 1.
    TDS_APPROX = 0
    Z_s = 0.0
    J_dip = array([1.0, 0.0, 0.], 'D')
    r_dip = array([0.1, 0.1, 20.0], 'd')
    list_of_test_edges_numbers = arange(N_RWG,dtype='i')
    list_of_src_edges_numbers = arange(N_RWG,dtype='i')
    MOM_FULL_PRECISION = 1
    EXCITATION = 'dipole'
    #CHOICE = "CFIE testing"
    #CHOICE = "fields verification"
    CHOICE = "dielectric target"
    if CHOICE=="CFIE testing":
        for coeff in [1.0, .8, 0.5, 0.2, 0.0]:
        #for coeff in [0.2]:
            #CFIE = array([coeff, coeff, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            CFIE = array([coeff, 0, -(1.0 - coeff) * sqrt(mu_0/eps_0), -(1.0 - coeff) * sqrt(mu_0/eps_0)], 'D')
            print
            print "CFIE =", CFIE
            target_MoM = Target_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r, mu_r, TDS_APPROX, Z_s, MOM_FULL_PRECISION)
            # excitation computation
            target_MoM.V_EH_computation(CFIE, target_mesh, J_dip, r_dip, w, eps_r, mu_r, list_of_test_edges_numbers, EXCITATION)
            #target_MoM.solveByInversion()
            target_MoM.solveByLUdecomposition()
            print "inverted MoM RCS =", sum(target_MoM.I_CFIE*target_MoM.V_EH[:,0])
            #computeCurrentsVisualization(w, target_mesh, target_MoM.I_CFIE)
            # now we try the iterative method
            count = 0
            I_CFIE_bicgstab = bicgstab(target_MoM.Z_CFIE_J, target_MoM.V_CFIE, x0=None, tol=1.0e-05, maxiter=1000, xtype=None, callback=itercount)
            print "bicgstab MoM RCS =", sum(I_CFIE_bicgstab[0]*target_MoM.V_EH[:,0]), "# of iterations =", count

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
        eps_r_in, mu_r_in = 3.0 - 0.0j, 1.0 - 0.0j
        k_out = w/c
        IS_PMCHWT = 1 # 1 if we want PMCHWT, 0 for CFIE
        TE, NE, TH, NH = 1., 1., 1., 1.
        CFIE = array([TE * coeff, NE * coeff, -TH * (1.0 - coeff) * sqrt(mu_0/(eps_0*eps_r_in)), -NH * (1.0 - coeff) * sqrt(mu_0/(eps_0*eps_r_in))], 'D')
        target_MoM = dielectricTarget_MoM(CFIE, list_of_test_edges_numbers, list_of_src_edges_numbers, target_mesh, w, eps_r_out, mu_r_out, eps_r_in, mu_r_in, MOM_FULL_PRECISION, IS_PMCHWT)
        N_J, N_M = len(list_of_test_edges_numbers), len(list_of_test_edges_numbers)
        
        
        #J_dip = array([0.0, 1.0, 0.], 'D') * sqrt((c/f)/sqrt(mu_0/eps_0) * 3.0/pi)
        #r_dip = array([0.0, 0.0, 0.2714], 'd')
        #target_MoM.V_EH_computation(CFIE, target_mesh, J_dip, r_dip, w, eps_r_out, mu_r_out, list_of_test_edges_numbers, EXCITATION)
        E_0 = array([-1.0,0.0,0.0],'D')
        k_hat = array([0.,0.,-1.0],'d')
        r_ref = zeros(3, 'd') #sum(triangles_centroids, axis=0)/T
        target_MoM.V_EH_plane_computation(CFIE, target_mesh, E_0, k_hat, r_ref, w, eps_r, mu_r, list_of_test_edges_numbers, IS_PMCHWT)
        #target_MoM.solveByInversion()
        target_MoM.solveByLUdecomposition()

        count = 0
        I_bicgstab = bicgstab(target_MoM.Z,target_MoM.V,x0=None,tol=1.0e-04, maxiter=1000,xtype=None, callback=itercount)[0]
        print "bicgstab MoM  # of iterations =", count

        count = 0
        I_gmres = gmres(target_MoM.Z, target_MoM.V, x0=None, tol=1.0e-04, maxiter=1000, xtype=None, restrt=300, callback=itercount)[0]
        print "gmres MoM  # of iterations =", count
        
        ## computation of the scattered fields
        J_obs = array([1.0, 0.0, 0.0], 'D')
        E_x, E_inc_x, E_tot_x = [], [], []
        for z in arange(-1.,-0.075,0.005):
            r_obs = array([0.0, 0.0, z], 'd')
            V_EH = computeV_EH(target_mesh, J_obs, r_obs, w, eps_r_out, mu_r_out, list_of_test_edges_numbers, 'dipole', 'F')
            E_x.append(sum(-target_MoM.I[:N_J]*V_EH[:,0] + target_MoM.I[N_J:N_J + N_M]*V_EH[:,2]))
            E_inc_x.append(E_0[0] * exp(-1.j * k_out * dot(k_hat, r_obs - r_ref)))
            E_tot_x.append(E_x[-1] + E_inc_x[-1])
            #print "Ex scatt inv =", E_z[-1], "z =", z 
            #print "Ex scatt gmres =", sum(-I_CFIE_gmres[:N_J]*V_EH[:,0] + I_CFIE_gmres[N_J:N_J + N_M]*V_EH[:,2]), "r_obs =", r_obs 
            #print "Ex scatt bicgstab =", sum(-I_CFIE_bicgstab[:N_J]*V_EH[:,0] + I_CFIE_bicgstab[N_J:N_J + N_M]*V_EH[:,2])
        for z in arange(-0.07,0.07,0.005):
            r_obs = array([0.0, 0.0, z], 'd')
            V_EH = computeV_EH(target_mesh, J_obs, r_obs, w, eps_r_in, mu_r_in, list_of_test_edges_numbers, 'dipole', 'F')
            E_x.append(sum(target_MoM.I[:N_J]*V_EH[:,0] - target_MoM.I[N_J:N_J + N_M]*V_EH[:,2]))
            E_inc_x.append(E_0[0] * exp(-1.j * k_out * dot(k_hat, r_obs - r_ref)))
            E_tot_x.append(E_x[-1])
        for z in arange(0.075,1.0,0.005):
            r_obs = array([0.0, 0.0, z], 'd')
            V_EH = computeV_EH(target_mesh, J_obs, r_obs, w, eps_r_out, mu_r_out, list_of_test_edges_numbers, 'dipole', 'F')
            E_x.append(sum(-target_MoM.I[:N_J]*V_EH[:,0] + target_MoM.I[N_J:N_J + N_M]*V_EH[:,2]))
            E_inc_x.append(E_0[0] * exp(-1.j * k_out * dot(k_hat, r_obs - r_ref)))
            E_tot_x.append(E_x[-1] + E_inc_x[-1])

        fr = open('E_tot_x_real.txt','w')
        fi = open('E_tot_x_imag.txt','w')
        for elem in E_tot_x:
            fr.write(str(real(elem)) + '\n')
            fi.write(str(imag(elem)) + '\n')
        fr.close()
        fi.close()
        from pylab import rc, plot, show, xlabel, ylabel, xticks, yticks, grid, legend, title
        #rc('text', usetex=True)
        FontSize=18
        LineWidth=2
        LEGEND = []
        plot(real(array(E_x)), 'bo-', linewidth = LineWidth)
        plot(real(array(E_tot_x)), 'go-', linewidth = LineWidth)
        LEGEND.append('E_x')
        LEGEND.append('E_tot_x')
        figureTitle = "Field E_x"
        figureTitle += ", f = " + str(f/1.e9) + " GHz"
        title(figureTitle,fontsize=FontSize+2)
        xlabel('z',fontsize=FontSize+2)
        ylabel('E_x',fontsize=FontSize+2)
        legend(LEGEND)
        xticks(fontsize=FontSize)
        yticks(fontsize=FontSize)
        grid(True)
        show()


        # incoming field computation
        #G_EJ, G_HJ = G_EJ_G_HJ(r_dip, r_obs, eps_r_out, mu_r_out, w/c)
        #E_inc, H_inc = dot(G_EJ, J_obs), dot(G_HJ, J_obs)
        #print "incoming E_field at observation point =", E_inc
        
        #print "Ex scatt inv =", sum(-target_MoM.I_CFIE[:N_J]*V_EH[:,0] + target_MoM.I_CFIE[N_J:N_J + N_M]*V_EH[:,2])
        #print "Ex scatt bicgstab =", sum(-I_CFIE_bicgstab[:N_J]*V_EH[:,0] + I_CFIE_bicgstab[N_J:N_J + N_M]*V_EH[:,2])
        #print "Ex scatt gmres =", sum(-I_CFIE_gmres[:N_J]*V_EH[:,0] + I_CFIE_gmres[N_J:N_J + N_M]*V_EH[:,2])
        #print "total =", sum(-target_MoM.I_CFIE[:N_J]*V_EH[:,0] + target_MoM.I_CFIE[N_J:N_J + N_M]*V_EH[:,2]) + E_inc[0]
        #print "bicgstab dielectric MoM RCS =", sum(-I_CFIE_iter[:N_J]*V_EH[:,0] + I_CFIE_iter[N_J:N_J + N_M]*V_EH[:,2]) + E_inc[0], "# of iterations =", target_MoM.iter_counter


