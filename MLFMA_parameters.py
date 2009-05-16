# import of the Params_Simu class
import sys, os
sys.path.append(os.path.abspath('./code'))
from Params_Simu import Params_Simu

# instanciation of Params_Simu class
params_simu = Params_Simu()

# do we use C++ or Python for mesh construction? Same but C++ is much more economic
languageForMeshConstruction = ["Python", "C++"]
params_simu.languageForMeshConstruction = languageForMeshConstruction[1]

##############################################################################
#
# Now the MoM/MLFMA parameters. Only for those who know what they're doing :)
#
##############################################################################
# COMPUTE_Z_NEAR = 1 or 0: do we need to compute the near field and SAI preconditioner 
# sparse matrices or not
params_simu.COMPUTE_Z_NEAR = 1

# MOM_FULL_PRECISION = 0/1: faster/slower Z_near computation but less/more precision
params_simu.MOM_FULL_PRECISION = 1

# the a (finest cubes sidelength) factor -- it will multiply lambda (the wavelength)
# to obtain the side length of the leaf (finest) level . Usually: a = lambda/4.
params_simu.a_factor = 0.25

# CFIE factor nu
params_simu.nu = 0.2

# BE_BH_N_Gauss_points: the number of points for the leaf cubes radiation functions calculation
# BE_BH_N_Gauss_points = 1, 3, 6, 9, 12, 13
params_simu.BE_BH_N_Gauss_points = 1

# do we allow parallelization by directions??
# this option should be enabled if there is a big number of processors wrt the number of cubes 
# at the ceiling level
params_simu.DIRECTIONS_PARALLELIZATION = 1

# do we allow a ceiling level that could be other than the third coarsest level?
# the ceiling level having the smallest memory footprint by its cubes is chosen
params_simu.ALLOW_CEILING_LEVEL = 0

# number of digits for the L computation
params_simu.NB_DIGITS = 3

# do we use the previous solution in monostatic computation?
params_simu.USE_PREVIOUS_SOLUTION = 1
# do we use the bistatic approximation for computing the monostatic RCS? 
# (much faster but less accurate if yes = 1)
params_simu.MONOSTATIC_BY_BISTATIC_APPROX = 0
params_simu.MAXIMUM_DELTA_PHASE = 0.0 # in degrees

# max block size for the near field and preconditioner matrices (in MBytes)
# the near field and preconditioner matrices are sliced in blocks and 
# dumped to the disk in order to minimize RAM memory occupation
params_simu.MAX_BLOCK_SIZE = 10.

# the integration type. Usually Gauss-Legendre for theta and Poncelet for phi.
params_simu.int_method_theta = "GAUSSL" 
params_simu.int_method_phi = "PONCELET"
params_simu.INCLUDE_BOUNDARIES = 0

# PERIODIC_XXX = 1 if we use Gaussian interpolator, 0 if we want to use the Lagrangian interpolator
# 1 is usually default here, as it yields better precision for the interpolator
params_simu.PERIODIC_Theta = 1
params_simu.PERIODIC_Phi = 1

# CYCLIC_XXX = 1 if we use cyclic interpolation for the radiation functions, 0 otherwise
# 1 is usually default here, as it yields better precision for the interpolator
params_simu.CYCLIC_Theta = 1
params_simu.CYCLIC_Phi = 1

# parameters for alpha translations smoothing and threshold
# this allows the sparsification of Alpha translation matrices
# Sparsification can be very important at coarse levels
# and thus we can have important memory savings for large targets
params_simu.alphaTranslation_smoothing_factor = 1.4
params_simu.alphaTranslation_thresholdRelValueMax = 1.0e-3
params_simu.alphaTranslation_RelativeCountAboveThreshold = 0.6

