# we import the base parameters
# if you want to modify them, you can do it here
# or in MLFMA_parameters.py
from MLFMA_parameters import *
from EM_constants import *

# the geometry or target we are going to use in the simulations.
# The list of targets can be viewed by looking up the geo directory
params_simu.pathToTarget = './geo'
params_simu.targetName = 'cubi'
# the dimensions of the target (if applicable)
params_simu.lx = .2
params_simu.ly = .2
params_simu.lz = .2
# the vertical offset of the target
params_simu.z_offset = 0.0
# the dimensions scaling factor (some targets are designed in other units than meters)
# if dimensions are in cm, params_simu.targetDimensions_scaling_factor = 0.01
# if dimensions are in mm, params_simu.targetDimensions_scaling_factor = 0.001
params_simu.targetDimensions_scaling_factor = 1.0

# frequency
params_simu.f = 2.12e9

# the lc (characteristic length) factor -- it will multiply lambda (the wavelength)
# to obtain the average edge length (lc) of the mesh. Usually: lc ~= lambda/10.
params_simu.lc_factor = 1.0/10.0

# CAD and meshing tools
# the following parameter tells if we want to mesh (True) on-the-fly with GMSH.
# If mesh comes from another program, like GIS, set it to False
params_simu.meshToMake = True
# the mesh format. Currently only 3 are supported.
MESH_FORMAT = ['GMSH', 'GIS', 'IDEAS']
# GMSH is the default meshing program
params_simu.meshFormat = MESH_FORMAT[0]
if params_simu.meshToMake:
    params_simu.meshFormat = MESH_FORMAT[0]

#relative epsilon and mu of the host medium
params_simu.eps_r = 1. + 0.j
params_simu.mu_r = 1. + 0.j

# computation settings. Bistatic computation is default.
# one can compute a bistatic excitation, followed by monostatic,
# followed by SAR imaging. Just be aware of the computation time.
params_simu.BISTATIC = 1
params_simu.MONOSTATIC_RCS = 0
params_simu.MONOSTATIC_SAR = 0
# the following is important only if params_simu.MONOSTATIC_RCS == 1
# or if params_simu.MONOSTATIC_SAR == 1
# can only be 1 (True) or 0 (False)
params_simu.COMPUTE_RCS_HH = 1
params_simu.COMPUTE_RCS_VV = 1
params_simu.COMPUTE_RCS_HV = 0
# safeguard against bad user choices
if (params_simu.COMPUTE_RCS_HH==0) and (params_simu.COMPUTE_RCS_VV==0):
    params_simu.COMPUTE_RCS_VV = 1

# Bistatic computation settings
# type, strenght and phase and direction, origin of the source.
# It will not be used in case of monostatic_RCS computation
EXCITATIONS = ['dipole', 'plane']
params_simu.EXCITATION = EXCITATIONS[0] # 0 is dipole, 1 is plane
if params_simu.EXCITATION=='dipole':
    # origin, strength, phase and polarization of the dipole
    params_simu.J_src_x = 1.0+0.j
    params_simu.J_src_y = 0.0+0.j
    params_simu.J_src_z = 0.0+0.j
    params_simu.r_src_x = 0.1
    params_simu.r_src_y = 0.1
    params_simu.r_src_z = 20.0
elif params_simu.EXCITATION=='plane':
    # origin, strength, phase and polarization of the plane wave
    params_simu.theta_inc = pi/2.0
    params_simu.phi_inc = pi/2.0
    params_simu.E_inc_theta = 1.0+0.j # the theta component
    params_simu.E_inc_phi = 0.0+0.j # the phi component
else:
    pass
# sampling points: sampling of the resulting field at user-specified points in space.
# It will be used only for BISTATIC
# The lists have to be of equal lengths.
# You can construct it by using a small program or list comprehension.
# example:
# params_simu.r_obs_x = [x1, x2, x3, ..., xn]
# params_simu.r_obs_y = [y1, y2, y3, ..., yn]
# params_simu.r_obs_z = [z1, z2, z3, ..., zn]
# will yield the points r1 = [x1, y1, z1], r2 = [x2, y2, z2], etc.
params_simu.r_obs_x = [0.1, 0.1]
params_simu.r_obs_y = [0.1, 0.1]
params_simu.r_obs_z = [20.0, 20.1]

# the angles for the monostatic RCS or the bistatic far-field data.
# Normally the code provides "best angles of observation", best
# from a "field spatial information for minimal sampling size" point of view.
# if you want the program to choose the best points, set AUTOMATIC to True.
# if you want to provide your own sampling points (because you need less/more angles,
# for example), set AUTOMATIC to False. If you chose 1 point only, and AUTOMATIC is False,
# then the point will correspond to START_THETA for thetas and START_PHI for phis.
# thetas
params_simu.START_THETA = pi/2.0
params_simu.STOP_THETA = pi
params_simu.AUTOMATIC_THETAS = False
params_simu.USER_DEFINED_NB_THETA = 1
# phis
params_simu.START_PHI = 0.0
params_simu.STOP_PHI = 2.0 * pi
params_simu.AUTOMATIC_PHIS = True
params_simu.USER_DEFINED_NB_PHI = 200
# now the monostatic SAR settings
if params_simu.MONOSTATIC_SAR==1:
    # dipole antenna will "fly" in a plane defined by local x and y axis
    # we also need an origin for the plane, and a distribution of observation points
    params_simu.SAR_local_x_hat = [-1.0, 0.0, 0.0]
    params_simu.SAR_local_y_hat = [0.0, 0.0, 1.0]
    params_simu.SAR_plane_origin = [0.0, 50.0, 0.0]
    # span of the scanning rectangle (in meters)
    params_simu.SAR_x_span = 200.
    params_simu.SAR_y_span = 200.
    # offset of the scanning rectangle wrt its origin
    params_simu.SAR_x_span_offset = -100.
    params_simu.SAR_y_span_offset = -100.
    # the number of points in each direction
    params_simu.SAR_N_x_points = 8
    params_simu.SAR_N_y_points = 3

# SOLVER DATA
# iterative solver tolerance
params_simu.TOL = 1.e-3
# iterative solver: BICGSTAB, GMRES, RGMRES or FGMRES
SOLVERS = ["BICGSTAB", "GMRES", "RGMRES", "FGMRES"]
params_simu.SOLVER = SOLVERS[0]
params_simu.MAXITER = 150
# RESTART is only for (F)GMRES
params_simu.RESTART = 30
# inner solver characteristics. will be used only if FGMRES is used
params_simu.INNER_SOLVER = SOLVERS[0]
params_simu.INNER_TOL = 0.25
params_simu.INNER_MAXITER = 15
params_simu.INNER_RESTART = 30
# preconditioner type. No choice here
params_simu.PRECOND = "FROB"

# figure that shows the far field or the RCS of the target
params_simu.SHOW_FIGURE = 0

# CURRENTS_VISUALIZATION tells if we want to create a view of the currents in
# GMSH or not. if 1, you can use GMSH to view the currents that have been created.
# open with GMSH the *.pos files that are in the ./geo directory
# only works if BISTATIC == 1
params_simu.CURRENTS_VISUALIZATION = 0

# now a parameter that tells if we want to be silent (little to no output)
params_simu.VERBOSE = 1

# do we have a thin dielectric sheet or impedance boundary condition?
params_simu.TDS_APPROX = 0 # 0 is default here
if params_simu.TDS_APPROX == 1:
    params_simu.nu = 1. # in this case we only use EFIE
    params_simu.SOLVER = "RGMRES" # BiCGSTAB does not work so good for this type of problem
    # refractive index N (from MEDGYESI-MITSCHANG and WANG, "Hybrid Solutions for Large-Impedance
    # Coated Bodies of Revolution", IEEE Trans. Ant. Prop., Vol 34, n. 11, November 1986)
    # see the article for "valid" values of N
    N = (10. - 1.j*0.2)
    # surface impedance Z_s
    params_simu.Z_s = 377.0 * 1.0/N
    # target name and dimensions
    params_simu.lx = 10.0 * c/(2.0 * pi * params_simu.f) # see MEDGYESI-MITSCHANG and WANG 86
else:
    params_simu.Z_s = 0.0 + 0.0j
