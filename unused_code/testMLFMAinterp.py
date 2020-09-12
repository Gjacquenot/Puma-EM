from MLFMA import *
from FMM import *
from MoM import *
from scipy.interpolate.interpolate import *
from interpolation import decimate_2D_lfi, decimate_2D_splint, gaussInterpMatrix
import copy
#from mayaviRender import renderVTK
from mayavi.tools import imv

def plotFtheta(F, theta, indPhi):
    """This function plots F(theta) for a given phi in two manners:
    1) complete solution from F(theta)
    2) for negative theta, use -F(|theta|, phi+pi)
    This allows us to see if equation (31) from "Optimal Interpolation of Radiated Fields over a Sphere", IEEE AP november 1991 is correct"""
    Ntheta, Nphi = F.shape[0], F.shape[1]
    indPhiPlusPi = indPhi + Nphi/2
    print(Nphi, indPhi, indPhiPlusPi)
    Falternative = zeros(Ntheta, Float)
    for i in range(Ntheta):
        if theta[i]<0.0:
            Falternative[i] = -F[Ntheta-i-1, indPhiPlusPi]
            print(theta[i], theta[Ntheta-i-1])
        else:
            Falternative[i] = F[i, indPhi]
    plot(theta, F[:, indPhi], 'bo-', theta, Falternative, 'rx-')
    grid(True)
    show()


def shiftExp(r1, r2, Xthetas, Xphis, k):
    """computes S = exp(-j*k*dot(k_hat,r2-r1)) for all theta and phi"""
    S = zeros((Xthetas.shape[0], Xphis.shape[0]), 'F')
    k_hat = zeros(3, Float)
    NTheta, NPhi = Xthetas.shape[0], Xphis.shape[0]
    for m in range(NTheta):
        cosTheta = cos(Xthetas[m])
        sinTheta = sin(Xthetas[m])
        for n in range(NPhi):
            k_hat = array([sinTheta*cos(Xphis[n]), sinTheta*sin(Xphis[n]), cosTheta])
            S[m,n] = exp( -1.j*k*sum(k_hat * (r2-r1)) )
    return S

def target_FMM_interpolation(target_FMM_L, target_FMM_TRUE, NOrderTheta, PERIODIC_Theta, CYCLIC_Theta, int_method_theta, NOrderPhi, PERIODIC_Phi, CYCLIC_Phi, int_method_phi):
    """construct a target_FMM by interpolating the values of in target_FMM_L"""
    target_FMM_interp = copy.deepcopy(target_FMM_TRUE)
    target_FMM_interp.B_tEJ_src[:,:,:,:] = 0.0
    target_FMM_interp.B_CFIE_test[:,:,:,:] = 0.0
    # interpolation
    for i in (target_FMM_interp.list_of_edges_numbers):
        for l in range(target_FMM_interp.B_tEJ_src.shape[1]):
            # src interp
            target_FMM_interp.B_tEJ_src[i,l] = decimate_2D_lfi(target_FMM_TRUE.Xtheta, target_FMM_L.Xtheta, 0., pi, NOrderTheta, PERIODIC_Theta, CYCLIC_Theta, target_FMM_TRUE.Xphi, target_FMM_L.Xphi, 0., 2.*pi, NOrderPhi, PERIODIC_Phi, CYCLIC_Phi, target_FMM_L.B_tEJ_src[i,l])
            #target_FMM_interp.B_tEJ_src[i,l] = gaussInterpMatrix(target_FMM_TRUE.Xtheta, 0., pi, target_FMM_L.Xtheta, NOrderTheta, CYCLIC_Theta, int_method_theta, target_FMM_L_1.Xphi, 0., 2.*pi, target_FMM_L.Xphi, NOrderPhi, CYCLIC_Phi, int_method_phi,  target_FMM_L.B_tEJ_src[i,l])
            target_FMM_interp.B_tEJ_src[i,l] *= shiftExp(target_FMM_L.target_mesh.edges_numbers_cubes_centroids[i], target_FMM_interp.target_mesh.edges_numbers_cubes_centroids[i], target_FMM_TRUE.Xtheta, target_FMM_TRUE.Xphi, target_FMM_interp.k)

            # test interp
            target_FMM_interp.B_CFIE_test[i,l] = decimate_2D_lfi(target_FMM_TRUE.Xtheta, target_FMM_L.Xtheta, 0., pi, NOrderTheta, PERIODIC_Theta, CYCLIC_Theta, target_FMM_TRUE.Xphi, target_FMM_L.Xphi, 0., 2.*pi, NOrderPhi, PERIODIC_Phi, CYCLIC_Phi, target_FMM_L.B_CFIE_test[i,l])
            #target_FMM_interp.B_CFIE_test[i,l] = gaussInterpMatrix(target_FMM_TRUE.Xtheta, 0., pi, target_FMM_L.Xtheta, NOrderTheta, CYCLIC_Theta, int_method_theta, target_FMM_TRUE.Xphi, 0., 2.*pi, target_FMM_L.Xphi, NOrderPhi, CYCLIC_Phi, int_method_phi,  target_FMM_L.B_CFIE_test[i,l])
            target_FMM_interp.B_CFIE_test[i,l] *= shiftExp(target_FMM_interp.target_mesh.edges_numbers_cubes_centroids[i], target_FMM_L.target_mesh.edges_numbers_cubes_centroids[i], target_FMM_TRUE.Xtheta, target_FMM_TRUE.Xphi, target_FMM_interp.k)

    if 1:
        gplt.surf(abs(target_FMM_TRUE.B_tEJ_src[0,0,:,:]-target_FMM_interp.B_tEJ_src[0,0,:,:]))
        #gplt.surf(abs(target_FMM_L_1.B_tEJ_src[0,0]))
        A = abs(target_FMM_TRUE.B_tEJ_src[0,0]-target_FMM_interp.B_tEJ_src[0,0])
        #A = abs(target_FMM_interp.B_tEJ_src[0,0])
        v = imv.surf(arange(target_FMM_TRUE.Xphi.shape[0]), arange(target_FMM_TRUE.Xtheta.shape[0]), A/max(max(A))*5)
    return target_FMM_interp


if __name__=="__main__":
    path = './geo'
    name = 'strip2.msh'
    z_offset = 0.0
    target_mesh = MeshClass(path, name, z_offset)

    CFIE = array([1, 0, 0, 0], Complex)
    f = 3.12e9
    w = 2. * pi * f
    eps_r = 1.
    mu_r = 1.
    k = w * sqrt(eps_0*eps_r*mu_0*mu_r) + 0.j
    a = 0.1*c/f
    J_dip = array([1., 0., 0.], Complex64)
    r_dip = array([0.1, 0.1, 1.], Float)

    target_mesh.cubes_data_computation(a)
    target_mesh.write_cubes(path, name)
    d0 = 6
    L  = L_computation(k, a)+4
    print("L =", L)
    N_points_theta = L+1
    N_points_phi = 2*L
    INCLUDE_BOUNDARIES = 0
    int_method_theta, int_method_phi = "GAUSSL", "TRAP"
    BE_BH_N_Gauss_points = 9
    MOM_FULL_PRECISION = 1 # for maximum/minimum precision but longer/shorter computation times

    list_of_edges_numbers = arange(max(target_mesh.edges_numbers)+1)
    I_PQ = ones(list_of_edges_numbers.shape[0], 'F')# + 1.j * rand(B_tEJ_src.shape[0])
    target_FMM_L = Target_FMM(CFIE, list_of_edges_numbers, L, target_mesh, w, eps_r, mu_r, N_points_theta, N_points_phi, int_method_theta, int_method_phi, INCLUDE_BOUNDARIES, J_dip, r_dip, 'F', MOM_FULL_PRECISION, BE_BH_N_Gauss_points)
    ZI_FMM_far_L = target_FMM_L.matvec_far(I_PQ)

    # we double the size of the cubes
    a *= 2
    L = L_computation(k, a)
    target_mesh.cubes_data_computation(a)
    target_mesh.write_cubes(path, name)
    N_points_theta, N_points_phi = L+1, 2*L
    d0 = 6
    print("L_FMM =", L)
    print("FMM: N_points_theta =", N_points_theta, ", N_points_phi =", N_points_phi)
    target_FMM_L_1 = Target_FMM(CFIE, list_of_edges_numbers, L, target_mesh, w, eps_r, mu_r, N_points_theta, N_points_phi, int_method_theta, int_method_phi, INCLUDE_BOUNDARIES, J_dip, r_dip, 'F', MOM_FULL_PRECISION, BE_BH_N_Gauss_points)
    ZI_FMM_far_L_1 = target_FMM_L_1.matvec_far(I_PQ)

    NOrderTheta, NOrderPhi = 5, 5
    CYCLIC_Theta, CYCLIC_Phi = 1, 1
    PERIODIC_Theta, PERIODIC_Phi = 1, 1
    target_FMM_interp = target_FMM_interpolation(target_FMM_L, target_FMM_L_1, NOrderTheta, PERIODIC_Theta, CYCLIC_Theta, int_method_theta, NOrderPhi, PERIODIC_Phi, CYCLIC_Phi, int_method_phi)
    #print(target_FMM_interp.B_tEJ_src[0,0])
    print("ZI_FMM_far_L =", ZI_FMM_far_L)
    print("ZI_FMM_far_L_1 =", ZI_FMM_far_L_1)
    print("ZI_FMM_interp =", target_FMM_interp.matvec_far(I_PQ))
    #if 0:
        #PhiInd = 15
        #coordInd = 0
        #patchInd = 0
        #subplot(211)
        #plot(target_FMM_L_1.Xtheta, real(target_FMM_L_1.B_tEJ_src[patchInd,coordInd,:,PhiInd]), 'go-', target_FMM_interp.Xtheta, real(target_FMM_interp.B_tEJ_src[patchInd,coordInd,:,PhiInd]), 'rx-')
        #grid(True)
        #subplot(212)
        #plot(target_FMM_L_1.Xtheta, imag(target_FMM_L_1.B_tEJ_src[patchInd,coordInd,:,PhiInd]), 'go-', target_FMM_L_1.Xtheta, imag(target_FMM_interp.B_tEJ_src[patchInd,coordInd,:,PhiInd]), 'rx-')
        ##plot(target_FMM_L.Xtheta, imag(target_FMM_L.B_tEJ_src[patchInd,coordInd,:,PhiInd/2]), 'bo-')
        #grid(True)
        #show()
