from mpi4py import *
#import commands
#import time, copy, pickle
from scipy import log, log10, ceil, floor, arccos, product, where, real, sqrt
from EM_constants import *
#from scipy.linalg import gmres, bicgstab
#from pylab import plot, rc, subplot, xlabel, ylabel, legend, xticks, yticks, grid, gca, setp, show, title
#from meshClass import MeshClass
#from V_EH import V_EH_dipole, V_EH_plane
#from integration import *
#from MoMPostProcessing import *
#from string import lower
#from PyGmsh import executeGmsh, write_geo
#from FMM_matvecs import matvec_FMM2
#from FMM_translation import *
#from FMM_Znear import#, Z_near_size_computation
#from FMM_precond import MgPreconditionerComputation
#from Z_MoM import Z_MoM

def L_computation(k, a, NB_DIGITS):
    """this function computes the number of expansion poles given the wavenumber k and the sidelength a"""
    L_tmp1 = a*real(k)*sqrt(3) + NB_DIGITS * (a*real(k)*sqrt(3))**(1./3.)
    L_tmp2 = a*real(k)*sqrt(3) + NB_DIGITS * log10(a*real(k)*sqrt(3) + pi)
    L_tmp3 = a*real(k)*sqrt(3) + 1 * log(a*real(k)*sqrt(3) + pi)
    L_tmp4 = a*real(k)*sqrt(3) + 1.8 * NB_DIGITS**(2.0/3.0) * (a*real(k)*sqrt(3))**(1.0/3.0)
    #my_id = MPI.COMM_WORLD.Get_rank()
    #if (my_id==0):
        #print "L_tmp1 =", L_tmp1, ", L_tmp2 =", L_tmp2, ", L_tmp3 =", L_tmp3, ", L_tmp4 =", L_tmp4
    L2 = where(L_tmp2<5., 5, int(ceil( L_tmp2 )))
    #return L2
    return int(floor( L_tmp2 ))
    #return int(ceil( L_tmp2 ))

#def Z_near(CFIE, target_mesh, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, MOM_FULL_PRECISION):
    #"""this function computes the near non-diagonal part of the MoM method"""
    #pathToSaveTo = "./tmp/Z_tmp/chunk0"
    #os.mkdir(pathToSaveTo)
    #C = target_mesh.cubes_centroids.shape[0]
    #Z_CFIE_near = zeros(target_mesh.N_near, ELEM_TYPE)
    #pq_array = zeros((target_mesh.N_near, 2), 'i')
    #startIndex = 0
    #for i in range(C):
        #Z_CFIE_near[startIndex:startIndex + target_mesh.N_nearPerCube[i]], pq_array[startIndex:startIndex + target_mesh.N_nearPerCube[i], :] = Z_nearPerCube(pathToSaveTo, CFIE, i, target_mesh, w, eps_r, mu_r, ELEM_TYPE, Z_TMP_ELEM_TYPE, MOM_FULL_PRECISION)
        #startIndex += target_mesh.N_nearPerCube[i]
        #sys.stdout.write("\r" + "Z_near computation. Percentage = %.4s" %str(i*100./C))
        #sys.stdout.flush()
    #print
    #return Z_CFIE_near, pq_array

#def B_EJ_B_HJ_computation(CFIE, list_of_edges_numbers, vertexes_coord, triangles_vertexes, triangles_edges_kinds, triangles_edges_numbers, triangles_edges_signs, triangles_edges_lengths, triangles_edges_opp_vertexes, triangles_normals, triangles_areas, edges_numbers_triangles, edges_numbers_cubes_centroids, IS_SRC, k, w, mu_r, XcosTheta, Xphi, N_Gauss_points, ELEM_TYPE):

    ### we first find the indexes of the triangles that will participate to the computations 
    #indexes_triangles = edges_numbers_triangles_indexes(list_of_edges_numbers, edges_numbers_triangles)

    ### creation of the local B vectors 
    #N_edges = list_of_edges_numbers.shape[0]
    #N_points_theta, N_points_phi = XcosTheta.shape[0], Xphi.shape[0]
    #B_EH = zeros((N_edges, 2, N_points_theta, N_points_phi), ELEM_TYPE)

    ### Now we must give "local" numbers to the edges, such that:
    ##
    ##     local_edge_number = edges_numbers_local_edges_numbers[edge_number]
    ##
    ##  This "edges_numbers_local_edges_numbers" must be passed to the C++ code.
    ##
    ##  The size of "edges_numbers_local_edges_numbers" is equal to the maximum edge number
    ##  encountered in the "triangles_edges_numbers[indexes_triangles]" 2-D array
    #list_of_encountered_edges_numbers = reshape( take(triangles_edges_numbers, indexes_triangles, axis=0), (1,-1) )[0,:]
    #edges_numbers_local_edges_numbers = ones(max(list_of_encountered_edges_numbers) + 1, 'i') * -1
    #put(edges_numbers_local_edges_numbers, list_of_edges_numbers, arange(N_edges))

    ## Note that only the nonborder edges having a positive number
    ## in "edges_numbers_local_edges_numbers" will participate to the B_tEJ computations
    #wrapping_code = """IT_theta_IT_phi_B_EJ_B_HJ (B_EH, CFIE, indexes_triangles, edges_numbers_local_edges_numbers, vertexes_coord, triangles_vertexes, triangles_edges_kinds, triangles_edges_numbers, triangles_edges_signs, triangles_edges_lengths, triangles_edges_opp_vertexes, triangles_normals, triangles_areas, edges_numbers_cubes_centroids, IS_SRC, k, w, mu_r, XcosTheta, Xphi, N_Gauss_points);
    #"""
    #weave.inline(wrapping_code,
                 #['B_EH', 'CFIE', 'indexes_triangles', 'edges_numbers_local_edges_numbers', 'vertexes_coord', 'triangles_vertexes', 'triangles_edges_kinds', 'triangles_edges_numbers', 'triangles_edges_signs', 'triangles_edges_lengths', 'triangles_edges_opp_vertexes', 'triangles_normals', 'triangles_areas', 'edges_numbers_cubes_centroids', 'IS_SRC', 'k', 'w', 'mu_r', 'XcosTheta', 'Xphi','N_Gauss_points'],
                 #type_converters = converters.blitz,
                 #include_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 #library_dirs = ['./MoM/', './MoM/amos/zbesh/'],
                 #libraries = ['MoM', 'AMOS', 'g2c'],
                 #headers = ['<iostream>','<complex>','"FMM.h"'],
                 #compiler = 'gcc')
    #return B_EH

#def Weights2D_computation(Wtheta, Wphi, ELEM_TYPE):
    #N_points_theta, N_points_phi = Wtheta.shape[0], Wphi.shape[0]
    #Weights2D = ones((N_points_theta, N_points_phi), 'f')
    #for i in range(N_points_theta):
        #Weights2D[i] = Wtheta[i] * Wphi
    #return real(Weights2D.astype(ELEM_TYPE))


#def psolveDiagonal(self, b):
    #"""simple diagonal preconditioning"""
    #return self.diagPrecond * b

#def psolveLeftFrob(self, b):
    #"""Mg static pattern geometric preconditioner"""
    #return matvec_Z_PQ_near(self.M_LeftFrob, b, self.M_LeftFrob_pq_array, self.ELEM_TYPE)

#class Target_FMM:
    #def __init__(self, CFIE, list_of_edges_numbers, L, target_mesh, w, eps_r, mu_r, N_points_theta,  N_points_phi, int_method_theta, int_method_phi, INCLUDE_BOUNDARIES, J_dip, r_dip, ELEM_TYPE, Z_TMP_ELEM_TYPE, MOM_FULL_PRECISION, BE_BH_N_Gauss_points, EXCITATION):
        ## cleaning and re-creating ./tmp, a temp dir used for calculations
        #commands.getoutput("rm -rf " + os.path.abspath(os.path.join('.','tmp')))
        #os.mkdir( os.path.abspath(os.path.join('.','tmp')) )
        #os.mkdir( os.path.abspath(os.path.join('.','tmp/Z_tmp')) )
        #os.mkdir( os.path.abspath(os.path.join('.','tmp/Z_near')) )
        #os.mkdir( os.path.abspath(os.path.join('.','tmp/Mg_LeftFrob')) )
        ### the wavenumber
        #self.target_mesh = target_mesh
        #self.ELEM_TYPE = ELEM_TYPE
        #self.Z_TMP_ELEM_TYPE = Z_TMP_ELEM_TYPE
        #self.MOM_FULL_PRECISION = MOM_FULL_PRECISION
        #self.BE_BH_N_Gauss_points = BE_BH_N_Gauss_points
        #self.k = w * sqrt(eps_0*eps_r*mu_0*mu_r)
        #self.N_points_theta,  self.N_points_phi = N_points_theta, N_points_phi
        #self.int_method_theta, self.int_method_phi = int_method_theta, int_method_phi
        ## we take the abscissas for cos(theta) due to the better performance in the integration
        ## so in effect the variable of integration is cos(theta) and not theta itself
        #self.XcosTheta, self.Wtheta = integr_1D_X_W(-1., 1., N_points_theta, int_method_theta, INCLUDE_BOUNDARIES)
        ## FMM.cpp takes an array of thetas, so we have to arccos() the cos(theta) array
        #self.Xtheta = arccos(self.XcosTheta)[-1::-1]
        #self.Xphi, self.Wphi = integr_1D_X_W(0., 2.*pi, N_points_phi, "PONCELET", INCLUDE_BOUNDARIES)
        #if int_method_phi=="TRAP":
            #deltaPhi = self.Xphi[1]-self.Xphi[0]
            #self.Xphi -= deltaPhi/2.
            #self.Xphi[0] = 0.
        #self.Xtheta, self.Wtheta, self.Xphi, self.Wphi = self.Xtheta.astype(lower(ELEM_TYPE)), self.Wtheta.astype(lower(ELEM_TYPE)), self.Xphi.astype(lower(ELEM_TYPE)), self.Wphi.astype(lower(ELEM_TYPE))
        #self.CFIE = CFIE
        ### the radiation functions of the RWG functions
        #self.B_tEJ_src = B_EJ_B_HJ_computation(self.CFIE, list_of_edges_numbers, target_mesh.vertexes_coord, target_mesh.triangles_vertexes, target_mesh.triangles_edges_kinds, target_mesh.triangles_edges_numbers, target_mesh.triangles_edges_signs, target_mesh.triangles_edges_lengths, target_mesh.triangles_edges_opp_vertexes, target_mesh.triangles_normals, target_mesh.triangles_areas, target_mesh.edges_numbers_triangles, target_mesh.edges_numbers_cubes_centroids, 1, self.k, w, mu_r, self.Xtheta, self.Xphi, self.BE_BH_N_Gauss_points, self.ELEM_TYPE)
        ### the receiving functions of the RWG functions
        #self.B_CFIE_test = B_EJ_B_HJ_computation(self.CFIE, list_of_edges_numbers, target_mesh.vertexes_coord, target_mesh.triangles_vertexes, target_mesh.triangles_edges_kinds, target_mesh.triangles_edges_numbers, target_mesh.triangles_edges_signs, target_mesh.triangles_edges_lengths, target_mesh.triangles_edges_opp_vertexes, target_mesh.triangles_normals, target_mesh.triangles_areas, target_mesh.edges_numbers_triangles, target_mesh.edges_numbers_cubes_centroids, 0, self.k, w, mu_r, self.Xtheta, self.Xphi, self.BE_BH_N_Gauss_points, self.ELEM_TYPE)
        #self.Weights2D = Weights2D_computation(self.Wtheta, self.Wphi, self.ELEM_TYPE)
        #sys.stdout.write("size of B_tE_J_src = %.5s" %str(prod(self.B_tEJ_src.shape) * self.B_tEJ_src.itemsize/(1024.0**2)) + " MBytes \n")
        #sys.stdout.write("size of B_CFIE_test = %.5s" %str(prod(self.B_CFIE_test.shape) * self.B_CFIE_test.itemsize/(1024.0**2)) + " MBytes \n")
        #self.alpha = alpha_computation(target_mesh.cubes_centroids, target_mesh.a, self.k, L, self.Xtheta, self.Xphi, self.ELEM_TYPE)
        #self.alpha_A = alpha_computation_alternative(target_mesh.cubes_centroids, target_mesh.a, self.k, L, self.Xtheta, self.Xphi, self.ELEM_TYPE)
        #print "alpha shape =", self.alpha.shape, ", alpha size =", prod(self.alpha.shape) *self.alpha.itemsize, "Bytes"
        #self.target_mesh.N_nearBlockDiag, self.target_mesh.N_near, self.target_mesh.N_nearPerCube = Z_near_size_computation(target_mesh.cubes_lists_edges_numbers, target_mesh.cubesNeighborsIndexes)
        #print "N_nearBlockDiag =", self.target_mesh.N_nearBlockDiag, ", N_near =", self.target_mesh.N_near
        #print "Z_CFIE_near computation."
        #self.Z_CFIE_near, self.pq_array = Z_near(CFIE, self.target_mesh, w, eps_r, mu_r, self.ELEM_TYPE, self.Z_TMP_ELEM_TYPE, self.MOM_FULL_PRECISION)
        #R_NORM_TYPE_1 = target_mesh.a
        #pathToReadFrom = "./tmp/Z_tmp"
        #self.M_LeftFrob, self.M_LeftFrob_pq_array = MgPreconditionerComputation(target_mesh, R_NORM_TYPE_1, ELEM_TYPE, Z_TMP_ELEM_TYPE, pathToReadFrom)
        #self.diagPrecond = zeros(self.B_CFIE_test.shape[0], complex)
        #for i in range(self.Z_CFIE_near.shape[0]):
            #if self.pq_array[i,0]==self.pq_array[i,1]:
                #self.diagPrecond[self.pq_array[i,0]] = 1.0/self.Z_CFIE_near[i]
	#sys.stdout.write("Size of Z_CFIE_near = %.5s" %str( prod(self.Z_CFIE_near.shape)*self.Z_CFIE_near.itemsize/(1024.0**2) ) + " MBytes \n")
	#sys.stdout.write("Z_CFIE_near pq_array size = %.5s" %str( prod(self.pq_array.shape)*self.pq_array.itemsize/(1024.0**2) ) + " MBytes \n")
	#sys.stdout.write("Size of M_LeftFrob = %.5s" %str( prod(self.M_LeftFrob.shape)*self.M_LeftFrob.itemsize/(1024.0**2) ) + " MBytes \n")
	#sys.stdout.write("M_LeftFrob pq_array size = %.5s" %str( prod(self.M_LeftFrob_pq_array.shape)*self.M_LeftFrob_pq_array.itemsize/(1024.0**2) ) + " MBytes \n")
        #self.numberOfIterations = 0
        ## dipole excitation
        #self.J_dip, self.r_dip = J_dip, r_dip
        #if EXCITATION=='dipole':
            #self.V_EH = V_EH_dipole(self.J_dip, self.r_dip, list_of_edges_numbers, target_mesh.edges_numbers_triangles, target_mesh.vertexes_coord, target_mesh.triangles_vertexes, target_mesh.triangles_edges_numbers, target_mesh.triangles_edges_kinds, target_mesh.triangles_edges_signs, target_mesh.triangles_edges_lengths, target_mesh.triangles_edges_opp_vertexes, target_mesh.triangles_normals, target_mesh.triangles_areas, w, eps_r, mu_r).astype(ELEM_TYPE)
        ## plane wave excitation: we try to have the same field as for dipole excitation
        #elif EXCITATION=='plane':
            #self.V_EH = V_EH_plane(self.J_dip, self.r_dip, list_of_edges_numbers, target_mesh.edges_numbers_triangles, target_mesh.vertexes_coord, target_mesh.triangles_vertexes, target_mesh.triangles_edges_numbers, target_mesh.triangles_edges_kinds, target_mesh.triangles_edges_signs, target_mesh.triangles_edges_lengths, target_mesh.triangles_edges_opp_vertexes, target_mesh.triangles_normals, target_mesh.triangles_areas, target_mesh.triangles_centroids, w, eps_r, mu_r).astype(ELEM_TYPE)
        #else:
            #print "ERROR: Wrong excitation setting. Exiting"
            #sys.exit(1)
        ## more little definitions
        #self.w, self.eps_r, self.mu_r, self.L = w, eps_r, mu_r, L
        #self.list_of_edges_numbers = list_of_edges_numbers
    ##psolve = psolveDiagonal
    #psolve = psolveLeftFrob
    ##matvec = matvec_FMM1
    #matvec = matvec_FMM2
    #matvec_far = matvec_far1
    #matvecNear = matvec_near

#if __name__=="__main__":
    ## now the meshing
    #path = './geo'
    #targetName = 'cube'
    #f = 2.12e9
    #write_geo(path, targetName, 'lc', c/f/10.0)
    #write_geo(path, targetName, 'lx', 0.03)
    #write_geo(path, targetName, 'ly', 0.03)
    #write_geo(path, targetName, 'lz', 0.03)
    #executeGmsh(path, targetName, 0)
    #z_offset = 0.0
    #target_mesh = MeshClass(path, targetName, z_offset)

    #J_dip = array([1., 0., 0.], complex)
    #r_dip = array([0.0, 50.0, 0.0], 'f')
    #EXCITATION = 'dipole'
    #w = 2. * pi * f
    #eps_r = 1.
    #mu_r = 1.
    #k = w * sqrt(eps_0*eps_r*mu_0*mu_r)
    #nu = 0.7
    #ELEM_TYPE = 'D' # we work with doubles. Pure floats don't work with linalg.gmres()
    #Z_TMP_ELEM_TYPE = 'D' # the type of the temporary Z arrays stored on the disk
    #MOM_FULL_PRECISION = 1 # for maximum/minimum precision but longer/shorter computation times
    #BE_BH_N_Gauss_points = 9
    #int_method_theta, int_method_phi = "GAUSSL", "TRAP"
    #INCLUDE_BOUNDARIES = 0
    #TOL = 1.0e-5 # gmres tolerance
    #RESTART = 20
    #IS_CLOSED_SURFACE = 1 * (target_mesh.E*2/3 == target_mesh.T)
    #print "closed surface =", IS_CLOSED_SURFACE
    #CFIE = array([nu, 0, 0, -1.j*(1.0-nu) * IS_CLOSED_SURFACE], complex)
    #a = 1.0 * 0.25*c/f
    #target_mesh.cubes_data_computation(a)
    #target_mesh.write_cubes()
    #NB_DIGITS = 1
    #L  = L_computation(k, a, NB_DIGITS)
    #print "L =", L

    ## number of points for the integration
    #N_points_theta = L+1
    #N_points_phi = 2*L
    #print "N_points_theta =", N_points_theta, ", N_points_phi =", N_points_phi

    #list_of_edges_numbers = arange(max(target_mesh.edges_numbers)+1)

    #t0 = time.time()
    #target_FMM = Target_FMM(CFIE, list_of_edges_numbers, L, target_mesh, w, eps_r, mu_r, N_points_theta, N_points_phi, int_method_theta, int_method_phi, INCLUDE_BOUNDARIES, J_dip, r_dip, ELEM_TYPE, Z_TMP_ELEM_TYPE, MOM_FULL_PRECISION, BE_BH_N_Gauss_points, EXCITATION)
    #V_CFIE = zeros(target_FMM.V_EH.shape[0], ELEM_TYPE)
    #for i in range(4):
        #V_CFIE += (target_FMM.V_EH[:, i] * CFIE[i]).astype(ELEM_TYPE)
    ##I_FMM = linalg.gmres(target_FMM, V_CFIE, restrt=RESTART, tol=TOL)[0]
    #I_FMM = linalg.bicgstab(target_FMM, V_CFIE, tol=TOL)[0]
    #print
    #print time.time() - t0, "seconds for FMM resolution"
    #print "FMM RCS =", sum(I_FMM*target_FMM.V_EH[:,0])
    #if 0 and (target_mesh.E<4e3):
        #print "computing MoM solution"
        #MOM_FULL_PRECISION = 1
        #t0 = time.time()
        #Z_CFIE = Z_MoM(CFIE, list_of_edges_numbers, list_of_edges_numbers, target_mesh.edges_numbers_triangles, target_mesh.vertexes_coord, target_mesh.triangles_vertexes, target_mesh.triangles_edges_numbers, target_mesh.triangles_edges_kinds, target_mesh.triangles_edges_signs, target_mesh.triangles_edges_lengths, target_mesh.triangles_edges_opp_vertexes, w, eps_r, mu_r, MOM_FULL_PRECISION)
        #I_MoM = dot(linalg.inv(Z_CFIE), V_CFIE)
        #print time.time() - t0, "seconds for MoM resolution"
        #print "MoM RCS =", sum(I_MoM*target_FMM.V_EH[:,0])
        #index_sortIreal = argsort(I_MoM.real)
        #index_sortIimag = argsort(I_MoM.imag)

        #rc('text', usetex=True)
        #FontSize=18
        #LineWidth=1
        #subplot(211)
        #plot(arange(I_MoM.shape[0]), take(I_MoM.real, index_sortIreal, axis=0), arange(I_MoM.shape[0]),  take(I_FMM.real, index_sortIreal, axis=0), 'r--', linewidth = LineWidth)
        #xlabel(r'Number of RWG function',fontsize=FontSize+2)
        #ylabel(r'$\mathfrac{Re}\left\{I\right\}$',fontsize=FontSize+2)
        #legend([r'MoM',r'FMM'])
        #xticks(fontsize=FontSize)
        #yticks(fontsize=FontSize)
        #grid(True)
        #subplot(212)
        #plot(arange(I_MoM.shape[0]), take(I_MoM.imag, index_sortIimag, axis=0), arange(I_MoM.shape[0]),  take(I_FMM.imag, index_sortIimag, axis=0), 'r--', linewidth = LineWidth)
        #xlabel(r'Number of RWG function',fontsize=FontSize+2)
        #ylabel(r'$\mathfrac{Im}\left\{I\right\}$',fontsize=FontSize+2)
        #legend([r'MoM',r'FMM'])
        #xticks(fontsize=FontSize)
        #yticks(fontsize=FontSize)
        #grid(True)
        #ltext = gca().get_legend().get_texts()
        #setp(ltext[0], fontsize = FontSize)
        #setp(ltext[1], fontsize = FontSize)
        #show()
