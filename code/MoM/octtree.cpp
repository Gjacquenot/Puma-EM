#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace std;

#include "octtree.h"
#include "readWriteBlitzArrayFromFile.h"
#include "interpolation.h"

/****************************************************************************/
/********************************** Octtree *********************************/
/****************************************************************************/

Octtree::Octtree(const string octtree_data_path, const blitz::Array<double, 2>& cubes_centroids, const int proc_id, const int num_procs)
{
  octtreeDataPath = octtree_data_path;
  // verbose or not?
  readIntFromASCIIFile(octtree_data_path + "VERBOSE.txt", VERBOSE);
  int CUBES_DISTRIBUTION;
  readIntFromASCIIFile(octtree_data_path + "CUBES_DISTRIBUTION.txt", CUBES_DISTRIBUTION);

  this->setProcNumber(proc_id);
  this->setTotalNumProcs(num_procs);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "creating the tree on process " << proc_id << " from disk data" << endl;
  numberOfUpdates = 0;
  readIntFromASCIIFile(octtree_data_path + "N_active_levels.txt", N_levels);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "N active levels = " << N_levels << endl;
  readIntFromASCIIFile(octtree_data_path + "ALLOW_CEILING_LEVEL.txt", ALLOW_CEILING_LEVEL);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "ALLOW_CEILING_LEVEL = " << ALLOW_CEILING_LEVEL << endl;
  readIntFromASCIIFile(octtree_data_path + "DIRECTIONS_PARALLELIZATION.txt", DIRECTIONS_PARALLELIZATION);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "DIRECTIONS_PARALLELIZATION = " << DIRECTIONS_PARALLELIZATION << endl;
  readIntFromASCIIFile(octtree_data_path + "N_GaussOnTriangle.txt", N_GaussOnTriangle);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "N_GaussOnTriangle = " << N_GaussOnTriangle << endl;

  // theta data
  float A_theta, B_theta;
  readFloatFromASCIIFile(octtree_data_path + "A_theta.txt", A_theta);
  readFloatFromASCIIFile(octtree_data_path + "B_theta.txt", B_theta);
  int INCLUDED_THETA_BOUNDARIES, PERIODIC_Theta, CYCLIC_Theta, NOrderInterpTheta;
  readIntFromASCIIFile(octtree_data_path + "INCLUDED_THETA_BOUNDARIES.txt", INCLUDED_THETA_BOUNDARIES);
  readIntFromASCIIFile(octtree_data_path + "PERIODIC_Theta.txt", PERIODIC_Theta);
  readIntFromASCIIFile(octtree_data_path + "CYCLIC_Theta.txt", CYCLIC_Theta);
  readIntFromASCIIFile(octtree_data_path + "NOrderInterpTheta.txt", NOrderInterpTheta);
  blitz::Array<int, 1> octtreeNthetas;
  readIntBlitzArray1DFromASCIIFile(octtree_data_path + "octtreeNthetas.txt", octtreeNthetas);
  blitz::Array<float, 2> octtreeXthetas, octtreeWthetas;
  readFloatBlitzArray2DFromASCIIFile(octtree_data_path + "octtreeXthetas.txt", octtreeXthetas);
  readFloatBlitzArray2DFromASCIIFile(octtree_data_path + "octtreeWthetas.txt", octtreeWthetas);

  // phi data
  float A_phi, B_phi;
  readFloatFromASCIIFile(octtree_data_path + "A_phi.txt", A_phi);
  readFloatFromASCIIFile(octtree_data_path + "B_phi.txt", B_phi);
  int INCLUDED_PHI_BOUNDARIES, PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi;
  readIntFromASCIIFile(octtree_data_path + "INCLUDED_PHI_BOUNDARIES.txt", INCLUDED_PHI_BOUNDARIES);
  readIntFromASCIIFile(octtree_data_path + "PERIODIC_Phi.txt", PERIODIC_Phi);
  readIntFromASCIIFile(octtree_data_path + "CYCLIC_Phi.txt", CYCLIC_Phi);
  readIntFromASCIIFile(octtree_data_path + "NOrderInterpPhi.txt", NOrderInterpPhi);
  blitz::Array<int, 1> octtreeNphis;
  readIntBlitzArray1DFromASCIIFile(octtree_data_path + "octtreeNphis.txt", octtreeNphis);
  blitz::Array<float, 2> octtreeXphis, octtreeWphis;
  readFloatBlitzArray2DFromASCIIFile(octtree_data_path + "octtreeXphis.txt", octtreeXphis);
  readFloatBlitzArray2DFromASCIIFile(octtree_data_path + "octtreeWphis.txt", octtreeWphis);

  const int N_coord = 2;

  double leaf_side_length;
  readDoubleFromASCIIFile(octtree_data_path + "leaf_side_length.txt", leaf_side_length);


  blitz::Array<int, 1> LExpansion;
  readIntBlitzArray1DFromASCIIFile(octtree_data_path + "LExpansion.txt", LExpansion);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "L expansions = " << LExpansion << endl;
  blitz::Array<std::complex<float>, 1> CFIEcoeffs;
  readComplexFloatBlitzArray1DFromASCIIFile(octtree_data_path + "CFIEcoeffs.txt", CFIEcoeffs);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "CFIE coeffs = " << CFIEcoeffs << endl;
  blitz::Array<double, 1> bigCubeLowerCoord(3), bigCubeCenterCoord(3);
  readDoubleBlitzArray1DFromASCIIFile(octtree_data_path + "big_cube_lower_coord.txt", bigCubeLowerCoord);
  readDoubleBlitzArray1DFromASCIIFile(octtree_data_path + "big_cube_center_coord.txt", bigCubeCenterCoord);
  for (int i=0 ; i<3 ; ++i) {
    big_cube_lower_coord[i] = bigCubeLowerCoord(i);
    big_cube_center_coord[i] = bigCubeCenterCoord(i);
  }


  readFloatFromASCIIFile(octtree_data_path + "w.txt", this->w);
  readComplexDoubleFromASCIIFile(octtree_data_path + "k.txt", this->k);
  readComplexFloatFromASCIIFile(octtree_data_path + "eps_r.txt", this->eps_r);
  readComplexFloatFromASCIIFile(octtree_data_path + "mu_r.txt", this->mu_r);
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "w = " << w << ", mu_rel = " << mu_r << ", eps_rel = " << eps_r << endl;
  CFIE.resize(CFIEcoeffs.size());
  CFIE = CFIEcoeffs;
  levels.reserve(N_levels); // reserve the memory
  float minCubesArraysSize;
  int indMinCubesArraysSize;
  { // construction of the first level
    blitz::Array<float, 1> XthetasNextLevel, XphisNextLevel;
    if (octtreeXthetas.extent(0) == 1) {
      XthetasNextLevel.resize(octtreeNthetas(0));
      XphisNextLevel.resize(octtreeNphis(0));
      XthetasNextLevel = octtreeXthetas(0, blitz::Range(0, octtreeNthetas(0)-1));
      XphisNextLevel = octtreeXphis(0, blitz::Range(0, octtreeNphis(0)-1));
    }
    else {
      XthetasNextLevel.resize(octtreeNthetas(1));
      XphisNextLevel.resize(octtreeNphis(1));
      XthetasNextLevel = octtreeXthetas(1, blitz::Range(0, octtreeNthetas(1)-1));
      XphisNextLevel = octtreeXphis(1, blitz::Range(0, octtreeNphis(1)-1));
    }
    levels.push_back( Level(N_levels+1,
                            LExpansion(0),
                            leaf_side_length,
                            this->big_cube_lower_coord,
                            cubes_centroids,
                            k,
                            N_coord,
                            A_theta,
                            B_theta,
                            octtreeXthetas(0, blitz::Range(0, octtreeNthetas(0)-1)),
                            octtreeWthetas(0, blitz::Range(0, octtreeNthetas(0)-1)),
                            INCLUDED_THETA_BOUNDARIES,
                            octtreeNthetas(0),
                            PERIODIC_Theta,
                            CYCLIC_Theta,
                            NOrderInterpTheta,
                            A_phi,
                            B_phi,
                            octtreeXphis(0, blitz::Range(0, octtreeNphis(0)-1)),
                            octtreeWphis(0, blitz::Range(0, octtreeNphis(0)-1)),
                            INCLUDED_PHI_BOUNDARIES,
                            octtreeNphis(0),
                            PERIODIC_Phi,
                            CYCLIC_Phi,
                            NOrderInterpPhi,
                            XthetasNextLevel,
                            XphisNextLevel,
                            VERBOSE ) );
    levels[0].sortCubesByParents();
    minCubesArraysSize = levels[0].getCubesSizeMB();
    indMinCubesArraysSize = 0;
  } // end of first level construction...using braces
  /// construction of the coarser levels
  for (int j=1 ; j<N_levels ; ++j) {
    blitz::Array<float, 1> XthetasNextLevel, XphisNextLevel;
    if (j<N_levels-1) {
      XthetasNextLevel.resize(octtreeNthetas(j+1));
      XthetasNextLevel = octtreeXthetas(j+1, blitz::Range(0, octtreeNthetas(j+1)-1));
      XphisNextLevel.resize(octtreeNphis(j+1));
      XphisNextLevel = octtreeXphis(j+1, blitz::Range(0, octtreeNphis(j+1)-1));
    }
    else {
      XthetasNextLevel.resize(octtreeNthetas(j));
      XthetasNextLevel = octtreeXthetas(j, blitz::Range(0, octtreeNthetas(j)-1));
      XphisNextLevel.resize(octtreeNphis(j));
      XphisNextLevel = octtreeXphis(j, blitz::Range(0, octtreeNphis(j)-1));
    }
    levels.push_back( Level( levels[j-1],
                             LExpansion(j),
                             this->big_cube_lower_coord,
                             N_coord,
                             A_theta,
                             B_theta,
                             octtreeXthetas(j, blitz::Range(0, octtreeNthetas(j)-1)),
                             octtreeWthetas(j, blitz::Range(0, octtreeNthetas(j)-1)),
                             INCLUDED_THETA_BOUNDARIES,
                             octtreeNthetas(j),
                             PERIODIC_Theta,
                             CYCLIC_Theta,
                             NOrderInterpTheta,
                             A_phi,
                             B_phi,
                             octtreeXphis(j, blitz::Range(0, octtreeNphis(j)-1)),
                             octtreeWphis(j, blitz::Range(0, octtreeNphis(j)-1)),
                             INCLUDED_PHI_BOUNDARIES,
                             octtreeNphis(j),
                             PERIODIC_Phi,
                             CYCLIC_Phi,
                             NOrderInterpPhi,
                             XthetasNextLevel,
                             XphisNextLevel,
                             VERBOSE ) );
    levels[j].sortCubesByParents();
    levels[j-1].updateFatherIndexes(levels[j]);
    if (levels[j].getCubesSizeMB() < minCubesArraysSize) {
      minCubesArraysSize = levels[j].getCubesSizeMB();
      indMinCubesArraysSize = j;
    }
    if (j==N_levels-1) levels[j].setCeiling(1); // we have reached the top level
  }
  if (CUBES_DISTRIBUTION==1) {
    // we are distributing the leaf cubes among the processors
    // we want a fair distribution, i.e. leaf cubes evenly partitioned between processes
    // therefore, we have to have a top level having a sufficient number of cubes
    // assignation for the near field computation
    int Z_NEAR_CUBES_DISTRIBUTION = 1;
    while ( (levels[(levels.size() - 1)].getLevelSize() < 25*num_procs) && (levels.size()>1) ) levels.pop_back();
    N_levels = levels.size();
    levels[N_levels-1].setCeiling(1);
    if ( (proc_id==0) && (VERBOSE==1) ) cout << "level chosen for Z_NEAR leaf cubes distribution = " << levels[N_levels-1].getLevel() << endl;
    assignCubesToProcessors(num_procs, Z_NEAR_CUBES_DISTRIBUTION);
    writeAssignedLeafCubesToDisk(octtree_data_path, "cubesIndexAndNumberToProcessNumber_FOR_Z_NEAR.txt");
  }
  else if (ALLOW_CEILING_LEVEL==1) { //we remove the non-necessary levels
    int j = levels.size() - 1;
    while (j>indMinCubesArraysSize) {
      levels.pop_back();
      j--;
    }
    levels[indMinCubesArraysSize].setCeiling(1);
    std::vector<Level>(levels).swap(levels); // trick for trimming exceeding capacity
    N_levels = levels.size();
  }
  if (CUBES_DISTRIBUTION==0) {
    for (int j=1 ; j<N_levels ; ++j) levels[j].searchCubesNeighborsIndexes();
    assignCubesToProcessors(num_procs, CUBES_DISTRIBUTION);
    for (int j=0 ; j<N_levels ; ++j) levels[j].computeLocalCubesIndexes(this->getProcNumber());
    // alpha translations 
    if ( (proc_id==0) && (VERBOSE==1) ) cout << "Searching the indexes of the possible cubes for alpha translations.........." << endl;
    for (int l=0; l<N_levels; ++l) findAlphaTransParticipantsIndexes(l);
    // level reduction computation...
    for (int j=0 ; j<N_levels ; ++j) levels[j].computeLevelReduction();
  }
  N_levels = levels.size();
  if ( (proc_id==0) && (VERBOSE==1) ) cout << "Tree construction terminated." << endl;
}

void Octtree::assignCubesToProcessors(const int num_procs, const int CUBES_DISTRIBUTION)
{
  // we use a top-down approach, which ensures that each cube
  // has all its descendants on the same processor
  int L = N_levels-1, NCubes, my_id = MPI::COMM_WORLD.Get_rank();
  if (N_levels==1) this->DIRECTIONS_PARALLELIZATION = 0; // get rid of a border case
  if ( (this->DIRECTIONS_PARALLELIZATION==1)&&(CUBES_DISTRIBUTION==0) ) { // we parallelize the last level by directions
    levels[L].DIRECTIONS_PARALLELIZATION = 1;
    NCubes = levels[L].getLevelSize();
    for (int i=0 ; i<NCubes ; ++i) levels[L].cubes[i].procNumber = -1; // -1 means cube belongs to all processes
    if ( (my_id==0) && (VERBOSE==1) ) {
      cout << "Process " << my_id << ". The total number of directions at level " << levels[L].getLevel() << " is N_directions = " << sum(levels[L].MPI_Scatterv_scounts) << endl;
      cout << "Process " << my_id << ". The sharing of directions between processes is as follows:" << endl; 
      cout << "Process " << my_id << "    levels[L].MPI_Scatterv_scounts = " << levels[L].MPI_Scatterv_scounts << endl;
      cout << "Process " << my_id << "    levels[L].MPI_Scatterv_displs = "<< levels[L].MPI_Scatterv_displs << endl;
    }
    L = L-1;
  }
  NCubes = levels[L].getLevelSize();
  if (NCubes<num_procs) {
    cout << "ERROR!! Octtree::assignCubesToProcessors: too few top-level cubes for the given number of processors" << endl;
    exit(1);
  }
  // top-most cell-parallelized level cubes attribution: bean-packing algorithm
  int procNumber = 0;
  for (int i=0 ; i<NCubes ; ++i) {
    levels[L].cubes[i].procNumber = procNumber;
    if (procNumber<num_procs-1) procNumber++; // we assign the next cube to the next process
    else procNumber = 0; // we start over again
  }
  // we need to set the sonsProcNumbers for the finest directions-parallelized level L+1
  // it is not important for the higher directions-parallelized levels
  if ( (this->DIRECTIONS_PARALLELIZATION==1)&&(CUBES_DISTRIBUTION==0) ) {
    NCubes = levels[L+1].getLevelSize();
    for (int i=0 ; i<NCubes ; ++i) {
      std::vector<int> sonsIndexes = levels[L+1].cubes[i].getSonsIndexes();
      const int N_sons = sonsIndexes.size();
      levels[L+1].cubes[i].sonsProcNumbers.resize(N_sons);
      for (int j=0 ; j<N_sons ; j++) {
        const int sonIndex = sonsIndexes[j];
        levels[L+1].cubes[i].sonsProcNumbers[j] = levels[L].cubes[sonIndex].procNumber;
      }
    }
  }
  // now we attribute the sons to processes for cell-parallelized levels
  // the sons have the same procNumber as their father
  for (int l=L ; l>0 ; l--) { // we go down the levels
    NCubes = levels[l].getLevelSize();
    for (int i=0 ; i<NCubes ; ++i) {
      std::vector<int> sonsIndexes = levels[l].cubes[i].getSonsIndexes();
      const int N_sons = sonsIndexes.size();
      levels[l].cubes[i].sonsProcNumbers.resize(N_sons);
      for (int j=0 ; j<N_sons ; j++) {
        const int sonIndex = sonsIndexes[j];
        levels[l-1].cubes[sonIndex].procNumber = levels[l].cubes[i].procNumber;
        levels[l-1].cubes[sonIndex].fatherProcNumber = levels[l].cubes[i].procNumber;
        levels[l].cubes[i].sonsProcNumbers[j] = levels[l-1].cubes[sonIndex].procNumber;
      }
    }
  }
}

void Octtree::writeAssignedLeafCubesToDisk(const string path, const string filename) {
  const int NC = levels[0].cubes.size();
  blitz::Array<int, 2> cubesNumberToProcessNumber(NC, 3);
  for (int i=0 ; i<NC ; ++i) {
    cubesNumberToProcessNumber(i, 0) = levels[0].cubes[i].getOldIndex();
    cubesNumberToProcessNumber(i, 1) = levels[0].cubes[i].getNumber();
    cubesNumberToProcessNumber(i, 2) = levels[0].cubes[i].getProcNumber();
  }
  writeIntBlitzArray2DToASCIIFile(path + filename, cubesNumberToProcessNumber);
}

void Octtree::computeIndexesOfCubesInOriginalMesh(blitz::Array<int, 1>& oldIndexesOfCubes) {
  levels[0].computeOldIndexesOfCubes(oldIndexesOfCubes);
}

void Octtree::computeGaussLocatedArguments(const blitz::Array<int, 1>& local_cubes_NRWG, 
                                           const blitz::Array<int, 1>& local_RWG_numbers, 
                                           const blitz::Array<int, 1>& local_RWG_Numbers_CFIE_OK, 
                                           const blitz::Array<int, 2>& local_RWGNumbers_signedTriangles, 
                                           const blitz::Array<float, 2>& local_RWGNumbers_trianglesCoord)
{
  if ( (getProcNumber()==0) && (VERBOSE==1) ) cout << "computing the leaf level Gauss located arguments.........." << endl;
  levels[0].computeGaussLocatedArguments(local_cubes_NRWG, local_RWG_numbers, local_RWG_Numbers_CFIE_OK, local_RWGNumbers_signedTriangles, local_RWGNumbers_trianglesCoord, this->N_GaussOnTriangle);
}

void Octtree::RWGs_renumbering(void)
{
  levels[0].RWGs_renumbering();
}

void Octtree::constructArrays(void) 
{
  const int N_levels = levels.size();
  if ( (getProcNumber()==0) && (VERBOSE==1) ) cout << "computing the alpha translations and shifting arrays.........." << endl;
  float alphaTranslation_smoothing_factor, alphaTranslation_thresholdRelValueMax, alphaTranslation_RelativeCountAboveThreshold;
  readFloatFromASCIIFile(this->octtreeDataPath + "alphaTranslation_smoothing_factor.txt", alphaTranslation_smoothing_factor);
  readFloatFromASCIIFile(this->octtreeDataPath + "alphaTranslation_thresholdRelValueMax.txt", alphaTranslation_thresholdRelValueMax);
  readFloatFromASCIIFile(this->octtreeDataPath + "alphaTranslation_RelativeCountAboveThreshold.txt", alphaTranslation_RelativeCountAboveThreshold);
  if (alphaTranslation_smoothing_factor<1.0) alphaTranslation_smoothing_factor = 1.0;
  if (alphaTranslation_smoothing_factor>2.0) alphaTranslation_smoothing_factor = 2.0;
  if (alphaTranslation_RelativeCountAboveThreshold > 1.0) alphaTranslation_RelativeCountAboveThreshold = 1.0;
  if (alphaTranslation_RelativeCountAboveThreshold < 0.0) alphaTranslation_RelativeCountAboveThreshold = 0.0;
  for (int j=0 ; j<N_levels ; ++j) {
    
    levels[j].NCubesXYZComputation(VERBOSE);
    levels[j].alphaTranslationsComputation(VERBOSE, alphaTranslation_smoothing_factor, alphaTranslation_thresholdRelValueMax, alphaTranslation_RelativeCountAboveThreshold);
    levels[j].shiftingArraysComputation();
  }
  if (getProcNumber()==0) {
    blitz::Array<double, 1> alphaTranslationsSizeMB(N_levels);
    blitz::Array<double, 1> cubesSizeMB(N_levels);
    blitz::Array<double, 1> shiftingArraysSizeMB(N_levels);
    for (int j=0 ; j<N_levels ; ++j) {
      alphaTranslationsSizeMB(j) = levels[j].getAlphaTranslationsSizeMB();
      cubesSizeMB(j) = levels[j].getCubesSizeMB();
      shiftingArraysSizeMB(j) = levels[j].getShiftingArraysSizeMB();
      if (VERBOSE==1) {
        cout << "for level " << levels[j].getLevel() << ", sizes of arrays are:" << endl;
        cout << "  size of alphaTranslations = " << alphaTranslationsSizeMB(j) << " MB;" << endl;
        cout << "  size of alphaTranslationsIndexes = " <<  levels[j].getAlphaTranslationsIndexesSizeMB() << " MB;" << endl;
        cout << "  size of cubes = " << cubesSizeMB(j) << " MB, number of cubes = " << levels[j].cubes.size() << endl;
        cout << "  size of indexes for alphaTranslations candidates = " << levels[j].getSizeMBOfAlphaTransParticipantsIndexes() << " MB;" << endl;
        cout << "  size of shiftingArrays = " << shiftingArraysSizeMB(j) << " MB;" << endl;
      }
    }
  }
}

void Octtree::copyOcttree(const Octtree& octtreeTocopy) /// copy constructor
{
  cout << "\ncopying octtree... " << endl;
  numberOfUpdates = octtreeTocopy.getNumberOfUpdates();
  procNumber = octtreeTocopy.getProcNumber();
  totalNumProcs = octtreeTocopy.getTotalNumProcs();
  L = octtreeTocopy.getL();
  k = octtreeTocopy.getK();
  N_GaussOnTriangle = octtreeTocopy.N_GaussOnTriangle;
  eps_r = octtreeTocopy.getEps_r();
  mu_r = octtreeTocopy.getMu_r();
  w = octtreeTocopy.getW();
  octtreeDataPath = octtreeTocopy.octtreeDataPath;
  int NLevels = octtreeTocopy.getLevelsSize();
  cout << "The number of Levels is " << NLevels << endl;
  levels.resize(NLevels);
  for (int j=0 ; j<NLevels ; j++) levels[j].copyLevel(octtreeTocopy.getLevel(j));
  for (int i=0 ; i<3 ; i++) big_cube_lower_coord[i] = octtreeTocopy.big_cube_lower_coord[i];
  for (int i=0 ; i<3 ; i++) big_cube_center_coord[i] = octtreeTocopy.big_cube_center_coord[i];
  CFIE.resize(octtreeTocopy.getCFIE().size());
  CFIE = octtreeTocopy.getCFIE();
  DIRECTIONS_PARALLELIZATION = octtreeTocopy.DIRECTIONS_PARALLELIZATION;
  N_levels = octtreeTocopy.N_levels;
  ALLOW_CEILING_LEVEL = octtreeTocopy.ALLOW_CEILING_LEVEL;
  VERBOSE = octtreeTocopy.VERBOSE;
  cout << "end of copying octtree... " << endl;
}

Octtree::Octtree(const Octtree& octtreeToCopy) /// copy constructor
{
  copyOcttree(octtreeToCopy);
}

Octtree& Octtree::operator=(const Octtree& octtreeToCopy) { /// copy assignment
  copyOcttree(octtreeToCopy);
  return *this;
}

Octtree::~Octtree()
{
  levels.clear();
  CFIE.free();
}

std::vector<int> Octtree::getNeighborsSonsIndexes(const int index, const int l) const
{
  std::vector<int> sonsOfNeighbors, sonsOfNeighborsTmp;
  const std::vector<int> neighborsIndexes(getCubeLevel(index, l).getNeighborsIndexes());
  for (unsigned int m=0 ; m<neighborsIndexes.size() ; ++m)
  {
    sonsOfNeighborsTmp = getCubeLevel(neighborsIndexes[m], l).getSonsIndexes();
    for (unsigned int n=0 ; n<sonsOfNeighborsTmp.size() ; ++n) sonsOfNeighbors.push_back(sonsOfNeighborsTmp[n]);
  }
  return sonsOfNeighbors;
}

void Octtree::findAlphaTransParticipantsIndexes(const int l)
/* we only fill in the cubes that are local, i.e. located on the same process as the octtree */
{
  //const int N_levels = levels.size();
  //for (int l=0; l<N_levels; ++l) {
    /* Each level has a "listOfFcToBeSent" and a "listOfFcToBeReceived".
     * Each list is as long as the number of processes,
     * and each element is accessed by the process number with which the communication will be 
     * established. Each element consists of a list of indexes corresponding to the radiation 
     * functions to be sent or received by the current process at the current level. 
     */
    std::vector< std::vector<int> > listOfFcToBeReceivedTmp, listOfFcToBeSentTmp;
    //levels[l].listOfFcToBeReceived.resize(this->getTotalNumProcs());
    //levels[l].listOfFcToBeSent.resize(this->getTotalNumProcs());
    listOfFcToBeReceivedTmp.resize(this->getTotalNumProcs());
    listOfFcToBeSentTmp.resize(this->getTotalNumProcs());
    
    std::vector<int> localCubesIndexes(levels[l].getLocalCubesIndexes());
    const int N_local_cubes = localCubesIndexes.size();
    const int N_cubes = levels[l].getLevelSize();
    int N_to_send(0), N_to_receive(0);
    // we treat the ceiling level differently than the regular levels
    // hereafter is the code for the ceiling level
    if (l==N_levels-1) { // if we are at the ceiling level
      for (int i=0; i<N_local_cubes; ++i) { // loop on the local cubes
        int indexLocalCube = localCubesIndexes[i];
        std::vector<int> localAlphaTransParticipantsIndexes, nonLocalAlphaTransParticipantsIndexes, nonLocalAlphaTransParticipantsProcNumbers;
        for (int j=0; j<N_cubes; ++j) { // loop on all the cubes (because ceiling level)
          const float * diffAbsCartCoord_1(levels[l].cubes[indexLocalCube].absoluteCartesianCoord); 
          const float * diffAbsCartCoord_2(levels[l].cubes[j].absoluteCartesianCoord);
          const float diffAbsCartCoord[3] = {diffAbsCartCoord_1[0]-diffAbsCartCoord_2[0], diffAbsCartCoord_1[1]-diffAbsCartCoord_2[1], diffAbsCartCoord_1[2]-diffAbsCartCoord_2[2]};
          bool condition = false;
          // condition = true if cubes are not touching
          for (int mm=0; mm<3; ++mm) condition = (condition || (abs(diffAbsCartCoord[mm]) > 1.0) );
          if (condition) {
            if ( levels[l].cubes[j].getProcNumber()==levels[l].cubes[indexLocalCube].getProcNumber() )
              localAlphaTransParticipantsIndexes.push_back(j);
            else { // the alphaTransPArticipant is not local
              nonLocalAlphaTransParticipantsIndexes.push_back(j);
              nonLocalAlphaTransParticipantsProcNumbers.push_back(j);
              listOfFcToBeReceivedTmp[levels[l].cubes[j].getProcNumber()].push_back(j);
              listOfFcToBeSentTmp[levels[l].cubes[j].getProcNumber()].push_back(indexLocalCube);
              N_to_send++;
              N_to_receive++;
            }
          }
        } // end for
        // we now trim the excess capacity of alphaTransParticipantsIndexes
        levels[l].cubes[indexLocalCube].localAlphaTransParticipantsIndexes.resize(localAlphaTransParticipantsIndexes.size());
        for (unsigned int j=0 ; j<localAlphaTransParticipantsIndexes.size() ; ++j) levels[l].cubes[indexLocalCube].localAlphaTransParticipantsIndexes[j] = localAlphaTransParticipantsIndexes[j];
        levels[l].cubes[indexLocalCube].nonLocalAlphaTransParticipantsIndexes.resize(nonLocalAlphaTransParticipantsIndexes.size());
        for (unsigned int j=0 ; j<nonLocalAlphaTransParticipantsIndexes.size() ; ++j) levels[l].cubes[indexLocalCube].nonLocalAlphaTransParticipantsIndexes[j] = nonLocalAlphaTransParticipantsIndexes[j];
      }
    }
    else { // for the NON-CEILING level
      for (int i=0; i<N_local_cubes; ++i) {
        int indexLocalCube = localCubesIndexes[i];
        std::vector<int> localAlphaTransParticipantsIndexes, nonLocalAlphaTransParticipantsIndexes, nonLocalAlphaTransParticipantsProcNumbers;
        const std::vector<int> possibleIndexes(getNeighborsSonsIndexes(levels[l].cubes[indexLocalCube].getFatherIndex(), l+1));
        for (unsigned int j=0; j<possibleIndexes.size(); j++) {// possible indexes of the alpha trans participants
          const int possibleIndex = possibleIndexes[j];
          const float * diffAbsCartCoord_1(levels[l].cubes[indexLocalCube].absoluteCartesianCoord);
          const float * diffAbsCartCoord_2(levels[l].cubes[possibleIndex].absoluteCartesianCoord);
          const float diffAbsCartCoord[3] = {diffAbsCartCoord_1[0]-diffAbsCartCoord_2[0], diffAbsCartCoord_1[1]-diffAbsCartCoord_2[1], diffAbsCartCoord_1[2]-diffAbsCartCoord_2[2]};
          bool condition = false;
          // condition = true if cubes are not touching
          for (int mm=0; mm<3; ++mm) condition = (condition || (abs(diffAbsCartCoord[mm]) > 1.0) );
          if (condition) {
            if ( levels[l].cubes[possibleIndex].getProcNumber()==levels[l].cubes[indexLocalCube].getProcNumber() )
              localAlphaTransParticipantsIndexes.push_back(possibleIndex);
            else { // the alphaTransPArticipant is not local
              nonLocalAlphaTransParticipantsIndexes.push_back(possibleIndex);
              nonLocalAlphaTransParticipantsProcNumbers.push_back(levels[l].cubes[possibleIndex].getProcNumber());
              listOfFcToBeReceivedTmp[levels[l].cubes[possibleIndex].getProcNumber()].push_back(possibleIndex);
              listOfFcToBeSentTmp[levels[l].cubes[possibleIndex].getProcNumber()].push_back(indexLocalCube);
              N_to_send++;
              N_to_receive++;
            }
          }
        }
        // we now trim the excess capacity of alphaTransParticipantsIndexes
        levels[l].cubes[indexLocalCube].localAlphaTransParticipantsIndexes.resize(localAlphaTransParticipantsIndexes.size());
        for (unsigned int j=0 ; j<localAlphaTransParticipantsIndexes.size() ; ++j) levels[l].cubes[indexLocalCube].localAlphaTransParticipantsIndexes[j] = localAlphaTransParticipantsIndexes[j];
        levels[l].cubes[indexLocalCube].nonLocalAlphaTransParticipantsIndexes.resize(nonLocalAlphaTransParticipantsIndexes.size());
        for (unsigned int j=0 ; j<nonLocalAlphaTransParticipantsIndexes.size() ; ++j) levels[l].cubes[indexLocalCube].nonLocalAlphaTransParticipantsIndexes[j] = nonLocalAlphaTransParticipantsIndexes[j];
      }
    }
    //if (VERBOSE==1) cout << "on process " << getProcNumber() << ", for level " << levels[l].getLevel() << ", the total number of radiation functions to be sent/received is N_to_send/N_to_receive = " << N_to_send << " / " << N_to_receive << endl;
    // we now eliminate the redundant radiations functions numbers in listOfFcToBeReceived
    // by sorting and eliminating redundant entries
    std::vector<int> listOfFcTmp;
    N_to_send = (N_to_receive = 0);
    for (int i=0 ; i<this->getTotalNumProcs() ; ++i) {
      if (i!= this->getProcNumber()) {
        if (listOfFcToBeReceivedTmp[i].size() > 0) {
          listOfFcTmp = listOfFcToBeReceivedTmp[i];
          sort(listOfFcTmp.begin(), listOfFcTmp.end());
          listOfFcToBeReceivedTmp[i].resize(0);
          listOfFcToBeReceivedTmp[i].push_back(listOfFcTmp[0]);
          for (unsigned int j=1 ; j<listOfFcTmp.size() ; ++j) {
            if (listOfFcTmp[j] != listOfFcToBeReceivedTmp[i].back()) listOfFcToBeReceivedTmp[i].push_back(listOfFcTmp[j]);
          }
          std::vector<int>(listOfFcToBeReceivedTmp[i]).swap(listOfFcToBeReceivedTmp[i]);
          N_to_receive += listOfFcToBeReceivedTmp[i].size();
        }
      }
    }
    levels[l].listOfFcToBeReceived = listOfFcToBeReceivedTmp;
    // we now eliminate the redundant radiations functions numbers in listOfFcToBeSent
    for (int i=0 ; i<listOfFcToBeSentTmp.size() ; ++i) {
      if (i!= this->getProcNumber()) {
        if (listOfFcToBeSentTmp[i].size() > 0) {
          listOfFcTmp = listOfFcToBeSentTmp[i];
          sort(listOfFcTmp.begin(), listOfFcTmp.end());
          listOfFcToBeSentTmp[i].resize(0);
          listOfFcToBeSentTmp[i].push_back(listOfFcTmp[0]);
          for (unsigned int j=1 ; j<listOfFcTmp.size() ; ++j) {
            if (listOfFcTmp[j] != listOfFcToBeSentTmp[i].back()) listOfFcToBeSentTmp[i].push_back(listOfFcTmp[j]);
          }
          std::vector<int>(listOfFcToBeSentTmp[i]).swap(listOfFcToBeSentTmp[i]);
          N_to_send += listOfFcToBeSentTmp[i].size();
        }
      }
    }
    levels[l].listOfFcToBeSent = listOfFcToBeSentTmp;
    if (VERBOSE==1) cout << "on process " << getProcNumber() << ", for level " << levels[l].getLevel() << ", the total number of Fc to be sent/received is N_to_send/N_to_receive = " << N_to_send << " / " << N_to_receive << endl;
  //} // for (int l=0; l<N_levels; ++l)
}

void Octtree::shiftExp(blitz::Array<std::complex<float>, 2> S,
                       const blitz::Array<std::complex<float>, 1>& shiftingArray)
{
  const int N = shiftingArray.size();
  //const int N_coord = 2;
  for (int i=0; i<N; i++) {
    S(0, i) *= shiftingArray(i);  
    S(1, i) *= shiftingArray(i);  
  }
}

void Octtree::S2DWeighting(blitz::Array<std::complex<float>, 2> S,
                           const blitz::Array<float, 1>& Wtheta,
                           const blitz::Array<float, 1>& Wphi)
{
  const int N_theta = Wtheta.size(), N_phi = Wphi.size();
  for (int i=0 ; i<N_theta ; ++i) {
    for (int j=0 ; j<N_phi ; ++j) {
      int index(i + j*N_theta);
      float W(Wtheta(i) * Wphi(j));
      S(0, index) *= W;
      S(1, index) *= W;
    }
  }
}

void Octtree::SupAlphaMultiplication(blitz::Array<std::complex<float>, 2>& SupAlpha,
                                     const blitz::Array<std::complex<float>, 2>& Sup,
                                     const blitz::Array<std::complex<float>, 1>& alphaTranslation,
                                     const blitz::Array<int, 1>& alphaTranslationIndexesNonZeros,
                                     const blitz::Array<int, 1>& alphaTranslationIndexes,
                                     const int alphaCartesianCoord[3])
{
  if ((abs(alphaCartesianCoord[0]) > 1) || (abs(alphaCartesianCoord[1]) > 1) || (abs(alphaCartesianCoord[2]) > 1)) {
    if (alphaTranslationIndexesNonZeros.size()==0) {
      const int N_alpha( alphaTranslationIndexes.size() );
      for (int i=0 ; i<N_alpha ; ++i) {
        const int newIndex(alphaTranslationIndexes(i));
        SupAlpha(0, i) += Sup(0, i) * alphaTranslation(newIndex);
        SupAlpha(1, i) += Sup(1, i) * alphaTranslation(newIndex);
      }
    }

    else {
      const int N_alpha(alphaTranslationIndexesNonZeros.size());
      for (int i=0 ; i<N_alpha ; ++i) {
        const int oldIndex = alphaTranslationIndexesNonZeros(i);
        const int newIndex = alphaTranslationIndexes(oldIndex);
        SupAlpha(0, newIndex) += Sup(0, newIndex) * alphaTranslation(i);
        SupAlpha(1, newIndex) += Sup(1, newIndex) * alphaTranslation(i);
      }
    }
  }
}

void Octtree::SupAlphaMultiplicationDirections(blitz::Array<std::complex<float>, 2>& SupAlpha,
                                               const blitz::Array<std::complex<float>, 2>& Sup,
                                               const blitz::Array<std::complex<float>, 1>& alphaTranslation,
                                               const blitz::Array<int, 1>& alphaTranslationIndexesNonZeros,
                                               const int alphaCartesianCoord[3])
{
  if ((abs(alphaCartesianCoord[0]) > 1) || (abs(alphaCartesianCoord[1]) > 1) || (abs(alphaCartesianCoord[2]) > 1)) {
    if ( (alphaTranslationIndexesNonZeros.size()==0) && (Sup.extent(1)==alphaTranslation.size()) ) {
      const int N_alpha = alphaTranslation.size();
      for (int i=0 ; i<N_alpha ; ++i) {
        SupAlpha(0, i) += Sup(0, i) * alphaTranslation(i);
        SupAlpha(1, i) += Sup(1, i) * alphaTranslation(i);
      }
    }
    else {
      const int N_alpha(alphaTranslationIndexesNonZeros.size());
      for (int i=0 ; i<N_alpha ; ++i) {
        const int index = alphaTranslationIndexesNonZeros(i);
        SupAlpha(0, index) += Sup(0, index) * alphaTranslation(i);
        SupAlpha(1, index) += Sup(1, index) * alphaTranslation(i);
      }
    }
  }
}

void Octtree::alphaTranslationsToCube(blitz::Array<std::complex<float>, 2>& S_tmp,
                                      const blitz::Array< blitz::Array<std::complex<float>, 2>, 1>& LevelSup,
                                      const int l,
                                      const int cubeIndex,
                                      const std::vector<int>& indexesAlphaParticipants,
                                      const int DIRECTIONS_PARALLELIZATION)
{
  blitz::Range all = blitz::Range::all();
  const float * cartCoord_1(levels[l].cubes[cubeIndex].absoluteCartesianCoord);
  S_tmp = 0.0;
  const int N_part = indexesAlphaParticipants.size();
  for (int j=0; j<N_part ; ++j) {
    const int indexParticipant = levels[l].cubesIndexesAfterReduction[indexesAlphaParticipants[j]];
    const float * cartCoord_2(levels[l].cubes[indexParticipant].absoluteCartesianCoord);
    const float DRcenters[3] = {cartCoord_1[0]-cartCoord_2[0], cartCoord_1[1]-cartCoord_2[1], cartCoord_1[2]-cartCoord_2[2]};
    const int alphaCartesianCoord[3] = { static_cast<int>( round(DRcenters[0]) ), static_cast<int>(  round(DRcenters[1]) ), static_cast<int>( round(DRcenters[2]) ) };
    if (DIRECTIONS_PARALLELIZATION!=1) {
      const int X = 1 * (alphaCartesianCoord[0]>=0), Y = 1 * (alphaCartesianCoord[1]>=0), Z = 1 * (alphaCartesianCoord[2]>=0);
      const int m = abs(alphaCartesianCoord[0]), n = abs(alphaCartesianCoord[1]), p = abs(alphaCartesianCoord[2]);
      SupAlphaMultiplication(S_tmp, LevelSup(indexParticipant), levels[l].alphaTranslations(m,n,p), levels[l].alphaTranslationsIndexesNonZeros(m, n, p), levels[l].alphaTranslationsIndexes(X, Y, Z, all), alphaCartesianCoord);
    }
    else {
      const int m = alphaCartesianCoord[0] + levels[l].getOffsetAlphaIndexX();
      const int n = alphaCartesianCoord[1] + levels[l].getOffsetAlphaIndexY();
      const int p = alphaCartesianCoord[2] + levels[l].getOffsetAlphaIndexZ();
      SupAlphaMultiplicationDirections(S_tmp, LevelSup(indexParticipant), levels[l].alphaTranslations(m,n,p), levels[l].alphaTranslationsIndexesNonZeros(m, n, p), alphaCartesianCoord);
    }
  }
}

void Octtree::updateSup(const blitz::Array<std::complex<float>, 1>& I_PQ) /// coefficients of RWG functions
{
  blitz::Range all = blitz::Range::all();
  int ierror;
  const int num_procs = getTotalNumProcs();
  numberOfUpdates += 1;
  if (this->getProcNumber()==0) cout << "\rTree update " << numberOfUpdates << ": level ";

  const int N_levels = levels.size(), my_id = this->getProcNumber();
  // first we update the Sups...which are stored (temporarily) in the Sdowns!!
  for (int l=0 ; l<N_levels ; ++l) {
    if (my_id==0) cout << levels[l].getLevel() << "..."; flush(cout);
    std::vector<int> localCubesIndexes(levels[l].getLocalCubesIndexes());
    const int N_local_cubes = localCubesIndexes.size();
    const int N_cubes = levels[l].getLevelSize();
    const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
    const int N_directions = (levels[l].DIRECTIONS_PARALLELIZATION!=1) ? N_theta*N_phi :  levels[l].MPI_Scatterv_scounts(my_id);
    if (levels[l].Sdown.size()==0) levels[l].Sdown.resize(N_cubes);
    // we first compute the Sups of all the cubes at the given level
    if (l==0) { // we use computeSup only for the leaf cubes
      for (int i=0 ; i<N_local_cubes ; ++i) {
        int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
        if (levels[l].Sdown(indexLocalCube).size()==0) levels[l].Sdown(indexLocalCube).resize(2, N_theta*N_phi);
        levels[l].computeSup(levels[l].Sdown(indexLocalCube), k, I_PQ, levels[l].cubes[indexLocalCube], levels[l].thetas, levels[l].phis);
      }
    }
    else if ( (l>0) && (levels[l].DIRECTIONS_PARALLELIZATION!=1) ) { // the Sups are obtained by interpolation from the sonsCubes Sups
      blitz::Array<std::complex<float>, 2> S_tmp(2, N_directions);
      for (int i=0 ; i<N_local_cubes ; ++i) {
        int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
        if (levels[l].Sdown(indexLocalCube).size()==0) levels[l].Sdown(indexLocalCube).resize(2, N_theta*N_phi);
        levels[l].Sdown(indexLocalCube) = 0.0;
        const int NSons = levels[l].cubes[indexLocalCube].sonsIndexes.size();
        for (int j=0 ; j<NSons ; ++j) {
          const int sonIndex = levels[l-1].cubesIndexesAfterReduction[levels[l].cubes[indexLocalCube].sonsIndexes[j]];
          // interpolation
          interpolate2Dlfi(S_tmp(0, all), levels[l-1].Sdown(sonIndex)(0, all), levels[l-1].lfi2D);
          interpolate2Dlfi(S_tmp(1, all), levels[l-1].Sdown(sonIndex)(1, all), levels[l-1].lfi2D);
          // shifting
          const float * rc_1(levels[l].cubes[indexLocalCube].rCenter);
          const float * rc_2(levels[l-1].cubes[sonIndex].rCenter);
          const double DRcenters[3] = {rc_1[0]-rc_2[0], rc_1[1]-rc_2[1], rc_1[2]-rc_2[2]};
          shiftExp( S_tmp, levels[l].getShiftingArray(DRcenters[0], DRcenters[1], DRcenters[2]) );
          levels[l].Sdown(indexLocalCube) += S_tmp;
        }
      }
    }
    else {
      const int sonLevel = l-1, N_directions = levels[l].MPI_Scatterv_scounts(my_id);
      std::vector<int> localCubesIndexesSonLevel(levels[sonLevel].getLocalCubesIndexes());
      int N_local_cubes_sonLevel = localCubesIndexesSonLevel.size(), max_N_local_cubes_sonLevel;
      blitz::Array<int, 1> FatherIndexes(num_procs), Array_N_local_cubes_sonLevel(num_procs);
      ierror = MPI_Allgather(&N_local_cubes_sonLevel, 1, MPI::INT, Array_N_local_cubes_sonLevel.data(), 1, MPI::INT, MPI_COMM_WORLD);
      max_N_local_cubes_sonLevel = max(Array_N_local_cubes_sonLevel);
      // initialization of the Sdowns Arrays
      for (int i=0 ; i<N_local_cubes ; ++i) {
        if (levels[l].Sdown(i).size()==0) levels[l].Sdown(i).resize(2, N_directions);
        levels[l].Sdown(i) = 0.0;
      }
      // now the aggregation stage
      blitz::Array<std::complex<float>, 2> S_tmp(2, N_theta*N_phi), S_tmp2(2, N_directions * num_procs);
      // we need to construct the receiving rcounts and rdispls arrays...
      blitz::Array<int, 1> scounts(levels[l].MPI_Scatterv_scounts), sdispls(levels[l].MPI_Scatterv_displs);
      blitz::Array<int, 1> rcounts(num_procs), rdispls(num_procs);
      for (int i=0 ; i<num_procs ; ++i) {
        rcounts(i) = scounts(my_id);
        rdispls(i) = i*scounts(my_id);
      }
      // now the serious loop...
      for (int i=0 ; i<max_N_local_cubes_sonLevel ; ++i) {
        int indexLocalCube = 0, fatherIndex = 0;
        if (i<N_local_cubes_sonLevel) {
          indexLocalCube = levels[sonLevel].cubesIndexesAfterReduction[localCubesIndexesSonLevel[i]];
          fatherIndex = levels[l].cubesIndexesAfterReduction[levels[sonLevel].cubes[indexLocalCube].getFatherIndex()];
          // interpolation
          interpolate2Dlfi(S_tmp(0, all), levels[sonLevel].Sdown(indexLocalCube)(0, all), levels[sonLevel].lfi2D);
          interpolate2Dlfi(S_tmp(1, all), levels[sonLevel].Sdown(indexLocalCube)(1, all), levels[sonLevel].lfi2D);
          // shifting
          const float * rc_1(levels[l].cubes[fatherIndex].rCenter);
          const float * rc_2(levels[sonLevel].cubes[indexLocalCube].rCenter);
          const double DRcenters[3] = {rc_1[0]-rc_2[0], rc_1[1]-rc_2[1], rc_1[2]-rc_2[2]};
          shiftExp( S_tmp, levels[l].getShiftingArray(DRcenters[0], DRcenters[1], DRcenters[2]) );
        }
        else S_tmp = 0.0;
        // we now "explode" the radiation functions...
        ierror = MPI_Alltoallv( S_tmp(0, all).data(), scounts.data(), sdispls.data(), MPI::COMPLEX, S_tmp2(0, all).data(), rcounts.data(), rdispls.data(), MPI::COMPLEX, MPI_COMM_WORLD);
        ierror = MPI_Alltoallv( S_tmp(1, all).data(), scounts.data(), sdispls.data(), MPI::COMPLEX, S_tmp2(1, all).data(), rcounts.data(), rdispls.data(), MPI::COMPLEX, MPI_COMM_WORLD);
        // we also need to pass the father indexes...
        ierror = MPI_Allgather(&fatherIndex, 1, MPI::INT, FatherIndexes.data(), 1, MPI::INT, MPI_COMM_WORLD);
        // we now perform aggregation at the parallelized-by-directions level...
        for (int j=0 ; j<FatherIndexes.size() ; ++j) {
          int fatherIndex = FatherIndexes(j);
          levels[l].Sdown(fatherIndex) += S_tmp2( all, blitz::Range(j*N_directions, (j+1)*N_directions - 1) );
        }
      }
    }
  }
  // now the translation stage!!!
  if (this->getProcNumber()==0) cout << "alpha Translations..."; flush(cout);
  for (int l=0 ; l<N_levels ; ++l) {
    std::vector<int> localCubesIndexes(levels[l].getLocalCubesIndexes());
    const int N_local_cubes = localCubesIndexes.size();
    const int N_cubes = levels[l].getLevelSize();
    const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
    const int N_directions = (levels[l].DIRECTIONS_PARALLELIZATION!=1) ? N_theta*N_phi :  levels[l].MPI_Scatterv_scounts(my_id);
    // we must define a SupThisLevel
    blitz::Array< blitz::Array<std::complex<float>, 2>, 1> SupThisLevel(N_cubes);
    for (int i=0 ; i<N_local_cubes ; ++i) {
      int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
      SupThisLevel(indexLocalCube).resize(2, N_directions);
      SupThisLevel(indexLocalCube) = levels[l].Sdown(indexLocalCube);
      levels[l].Sdown(indexLocalCube) = 0.0;
    }
    // alpha translations for cell-parallelized level
    if (levels[l].DIRECTIONS_PARALLELIZATION!=1) {
      if (l<3) exchangeSupsInBlocks(SupThisLevel, l, localCubesIndexes);
      else exchangeSupsIndividually(SupThisLevel, l, localCubesIndexes);
      // non local alpha translations
      blitz::Array<std::complex<float>, 2> S_tmp(2, N_theta*N_phi);
      for (int i=0 ; i<N_local_cubes ; ++i) {
        int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
        alphaTranslationsToCube(S_tmp, SupThisLevel, l, indexLocalCube, levels[l].cubes[indexLocalCube].getNonLocalAlphaTransParticipantsIndexes(), levels[l].DIRECTIONS_PARALLELIZATION);
        levels[l].Sdown(indexLocalCube) += S_tmp;
      }
    }
    // alpha translations for directions-parallelized level
    else {
      for (int i=0 ; i<N_local_cubes ; ++i) {
        int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
        alphaTranslationsToCube(levels[l].Sdown(indexLocalCube), SupThisLevel, l, indexLocalCube, levels[l].cubes[indexLocalCube].getLocalAlphaTransParticipantsIndexes(), levels[l].DIRECTIONS_PARALLELIZATION);
      }
    }
    ierror = MPI_Barrier(MPI::COMM_WORLD);
  }
}

void Octtree::exchangeSupsIndividually(blitz::Array< blitz::Array<std::complex<float>, 2>, 1>& SupThisLevel, const int l, const std::vector<int> & localCubesIndexes) {
  const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size(), N_coord = 2;
  std::vector< std::vector<MPI_Request> > isend_request, irecv_request;
  std::vector< std::vector<MPI_Status> > isend_status, irecv_status;
  isend_request.resize(getTotalNumProcs());
  irecv_request.resize(getTotalNumProcs());
  isend_status.resize(getTotalNumProcs());
  irecv_status.resize(getTotalNumProcs());
  int ierror = MPI_Barrier(MPI::COMM_WORLD);
  for (int i=0 ; i<getTotalNumProcs() ; ++i) {
    if (this->getProcNumber()!=i) {
      const int BUF_SIZE = N_theta * N_phi * N_coord;
      std::vector<int> listOfFcToBeSent(levels[l].getListOfFcToBeSent()[i]);
      std::vector<int> listOfFcToBeReceived(levels[l].getListOfFcToBeReceived()[i]);
      const int NToReceive = listOfFcToBeReceived.size();
      const int NToSend = listOfFcToBeSent.size();
      irecv_request[i].resize(NToReceive);
      isend_request[i].resize(NToSend);
      irecv_status[i].resize(NToReceive);
      isend_status[i].resize(NToSend);
      for (int j=0 ; j<NToReceive ; ++j) {
        const int indexRecv = levels[l].cubesIndexesAfterReduction[listOfFcToBeReceived[j]];
        SupThisLevel(indexRecv).resize(N_coord, N_theta*N_phi);
        ierror = MPI_Irecv(SupThisLevel(indexRecv).data(), BUF_SIZE, MPI::COMPLEX, i, listOfFcToBeReceived[j], MPI::COMM_WORLD, &irecv_request[i][j]);
      }
      for (int j=0 ; j<NToSend ; ++j) {
        const int indexSend = levels[l].cubesIndexesAfterReduction[listOfFcToBeSent[j]];
        ierror = MPI_Isend(SupThisLevel(indexSend).data(), BUF_SIZE, MPI::COMPLEX, i, listOfFcToBeSent[j], MPI::COMM_WORLD, &isend_request[i][j]);
      }
    }
  }
  // we then perform all the necessary local alpha translations at level l
  const int N_local_cubes = localCubesIndexes.size();
  for (int i=0 ; i<N_local_cubes ; ++i) {
    int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
    alphaTranslationsToCube(levels[l].Sdown(indexLocalCube), SupThisLevel, l, indexLocalCube, levels[l].cubes[indexLocalCube].getLocalAlphaTransParticipantsIndexes(), levels[l].DIRECTIONS_PARALLELIZATION);
  }
  // wait operation
  for (int i=0 ; i<getTotalNumProcs() ; ++i) {
    if (this->getProcNumber()!=i) {
      for (unsigned int j=0 ; j<irecv_request[i].size() ; ++j) ierror = MPI_Wait(&irecv_request[i][j], &irecv_status[i][j]);
      for (unsigned int j=0 ; j<isend_request[i].size() ; ++j) ierror = MPI_Wait(&isend_request[i][j], &isend_status[i][j]);
    }
  }
}

void Octtree::exchangeSupsInBlocks(blitz::Array< blitz::Array<std::complex<float>, 2>, 1>& SupThisLevel, const int l, const std::vector<int> & localCubesIndexes) {
  blitz::Range all = blitz::Range::all();
  const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size(), N_coord = 2;
  // we then communicate the necessary Fc radiation functions for alpha multiplication
  std::vector< MPI_Request > isend_request, irecv_request;
  std::vector< MPI_Status > isend_status, irecv_status;
  isend_request.resize(getTotalNumProcs());
  irecv_request.resize(getTotalNumProcs());
  isend_status.resize(getTotalNumProcs());
  irecv_status.resize(getTotalNumProcs());
  int ierror = MPI_Barrier(MPI::COMM_WORLD);
  // creation of the message buffers
  blitz::Array< blitz::Array<std::complex<float>, 2>, 1> buffToSend, buffToRecv;
  buffToSend.resize(getTotalNumProcs());
  buffToRecv.resize(getTotalNumProcs());
  for (int i=0 ; i<getTotalNumProcs() ; ++i) {
    buffToSend(i).resize(N_coord, levels[l].getListOfFcToBeSent()[i].size() * N_theta*N_phi);
    buffToRecv(i).resize(N_coord, levels[l].getListOfFcToBeReceived()[i].size() * N_theta*N_phi);
    // copying the Fc's to be sent in the send buffer
    std::vector<int> listOfFcToBeSent(levels[l].getListOfFcToBeSent()[i]);
    int startIndex = 0, stopIndex, FcIndex;
    for (int j=0 ; j<listOfFcToBeSent.size() ; ++j) {
      stopIndex = startIndex + N_theta*N_phi - 1;
      FcIndex = levels[l].cubesIndexesAfterReduction[listOfFcToBeSent[j]];
      buffToSend(i)(all, blitz::Range(startIndex, stopIndex)) = SupThisLevel(FcIndex);
      startIndex = stopIndex + 1;
    }
  }
  // and now sending and receiving
  for (int i=0 ; i<getTotalNumProcs() ; ++i) {
    if (this->getProcNumber()!=i) {
      int BUF_SIZE = buffToRecv(i).size(), FLAG = 22;
      ierror = MPI_Irecv(buffToRecv(i).data(), BUF_SIZE, MPI::COMPLEX, i, FLAG, MPI::COMM_WORLD, &irecv_request[i]);
    }
  }
  for (int i=0 ; i<getTotalNumProcs() ; ++i) {
    if (this->getProcNumber()!=i) {
      int BUF_SIZE = buffToSend(i).size(), FLAG = 22;
      ierror = MPI_Isend(buffToSend(i).data(), BUF_SIZE, MPI::COMPLEX, i, FLAG, MPI::COMM_WORLD, &isend_request[i]);
    }
  }
  // we then perform all the necessary alpha translations at level l
  const int N_local_cubes = localCubesIndexes.size();
  for (int i=0 ; i<N_local_cubes ; ++i) {
    int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
    alphaTranslationsToCube(levels[l].Sdown(indexLocalCube), SupThisLevel, l, indexLocalCube, levels[l].cubes[indexLocalCube].getLocalAlphaTransParticipantsIndexes(), levels[l].DIRECTIONS_PARALLELIZATION);
  }
  // wait operation
  for (int i=0 ; i<getTotalNumProcs() ; ++i) {
    if (this->getProcNumber()!=i) {
      ierror = MPI_Wait(&irecv_request[i], &irecv_status[i]);
      ierror = MPI_Wait(&isend_request[i], &isend_status[i]);
    }
  }
  buffToSend.free();
  // now we copy the received buffer to the Fc's
  for (int i=0 ; i<getTotalNumProcs() ; ++i) {
    if (this->getProcNumber()!=i) {
      std::vector<int> listOfFcToBeReceived(levels[l].getListOfFcToBeReceived()[i]);
      int startIndex = 0, stopIndex, FcIndex;
      for (int j=0 ; j<listOfFcToBeReceived.size() ; ++j) {
        stopIndex = startIndex + N_theta*N_phi - 1;
        FcIndex = levels[l].cubesIndexesAfterReduction[listOfFcToBeReceived[j]];
        SupThisLevel(FcIndex).resize(N_coord, N_theta * N_phi);
        SupThisLevel(FcIndex) = buffToRecv(i)(all, blitz::Range(startIndex, stopIndex));
        startIndex = stopIndex + 1;
      }
    }
  }
}

void Octtree::ZIFarComputation(blitz::Array<std::complex<float>, 1>& ZI, /// result of matrix-vector multiplication
                               const blitz::Array<std::complex<float>, 1>& I_PQ) /// coefficients of RWG functions
{
  // update of all the Sup of the tree
  blitz::Range all = blitz::Range::all();
  updateSup(I_PQ);
  if (this->getProcNumber()==0) cout << "tree descent"; flush(cout);
  const int N_levels = levels.size(), L = N_levels-1;
  int ierror, my_id = getProcNumber(), thisLevel = L;

  if (levels[thisLevel].DIRECTIONS_PARALLELIZATION==1) {
    const int sonLevel = thisLevel-1;
    std::vector<int> localCubesIndexes = levels[thisLevel].getLocalCubesIndexes();
    const int N_local_Cubes = localCubesIndexes.size();
    blitz::Array<std::complex<float>, 2> Stmp(2, sum(levels[thisLevel].MPI_Scatterv_scounts)), Stmp2(2, sum(levels[thisLevel].MPI_Scatterv_scounts)), Stmp3(2, levels[sonLevel].thetas.size() * levels[sonLevel].phis.size());
    for (int i=0 ; i<N_local_Cubes ; ++i) {
      int indexLocalCube = levels[thisLevel].cubesIndexesAfterReduction[localCubesIndexes[i]];
      ierror = MPI_Allgatherv ( levels[thisLevel].Sdown(indexLocalCube)(0, all).data(), levels[thisLevel].Sdown(indexLocalCube)(0, all).size(), MPI::COMPLEX, Stmp(0, all).data(), levels[thisLevel].MPI_Scatterv_scounts.data(), levels[thisLevel].MPI_Scatterv_displs.data(), MPI::COMPLEX, MPI::COMM_WORLD );
      ierror = MPI_Allgatherv ( levels[thisLevel].Sdown(indexLocalCube)(1, all).data(), levels[thisLevel].Sdown(indexLocalCube)(1, all).size(), MPI::COMPLEX, Stmp(1, all).data(), levels[thisLevel].MPI_Scatterv_scounts.data(), levels[thisLevel].MPI_Scatterv_displs.data(), MPI::COMPLEX, MPI::COMM_WORLD );

      std::vector<int> sonsIndexes = levels[thisLevel].cubes[indexLocalCube].sonsIndexes;
      std::vector<int> sonsProcNumbers = levels[thisLevel].cubes[indexLocalCube].sonsProcNumbers;
      for (unsigned int j=0 ; j<sonsIndexes.size() ; ++j) {
        // we need to do the following only if the sons are local
        if (my_id==sonsProcNumbers[j]) {
          const int sonIndex = levels[sonLevel].cubesIndexesAfterReduction[sonsIndexes[j]];
          // shifting
          const float * rc_1(levels[sonLevel].cubes[sonIndex].rCenter);
          const float * rc_2(levels[thisLevel].cubes[indexLocalCube].rCenter);
          const double DRcenters[3] = {rc_1[0]-rc_2[0], rc_1[1]-rc_2[1], rc_1[2]-rc_2[2]};
          Stmp2 = Stmp;
          shiftExp( Stmp2, levels[thisLevel].getShiftingArray(DRcenters[0], DRcenters[1], DRcenters[2]) );
          // anterpolate
          for (int m=0 ; m<2 ; m++) anterpolate2Dlfi(Stmp3(m, all), Stmp2(m, all), levels[sonLevel].lfi2D);
          levels[sonLevel].Sdown(sonIndex) += Stmp3;
        }
      }
    }
    thisLevel--;
  }

  // now the descent towards the bottom
  while (thisLevel>0) {
    const int sonLevel = thisLevel-1;
    std::vector<int> localCubesIndexes = levels[thisLevel].getLocalCubesIndexes();
    const int N_local_Cubes = localCubesIndexes.size();
    blitz::Array<std::complex<float>, 2> Stmp2(2, levels[thisLevel].thetas.size() * levels[thisLevel].phis.size()), Stmp3(2, levels[sonLevel].thetas.size() * levels[sonLevel].phis.size());
    for (int i=0 ; i<N_local_Cubes ; ++i) {
      int indexLocalCube = levels[thisLevel].cubesIndexesAfterReduction[localCubesIndexes[i]];
      std::vector<int> sonsIndexes = levels[thisLevel].cubes[indexLocalCube].sonsIndexes;
      for (unsigned int j=0 ; j<sonsIndexes.size() ; ++j) {
        const int sonIndex = levels[sonLevel].cubesIndexesAfterReduction[sonsIndexes[j]];
        // shifting
        const float * rc_1(levels[sonLevel].cubes[sonIndex].rCenter);
        const float * rc_2(levels[thisLevel].cubes[indexLocalCube].rCenter);
        const double DRcenters[3] = {rc_1[0]-rc_2[0], rc_1[1]-rc_2[1], rc_1[2]-rc_2[2]};
        Stmp2 = levels[thisLevel].Sdown(indexLocalCube);
        shiftExp( Stmp2, levels[thisLevel].getShiftingArray(DRcenters[0], DRcenters[1], DRcenters[2]) );
        // anterpolate
        for (int m=0 ; m<2 ; m++) anterpolate2Dlfi(Stmp3(m, all), Stmp2(m, all), levels[sonLevel].lfi2D);
        levels[sonLevel].Sdown(sonIndex) += Stmp3;
      }
    }
    thisLevel--;
  }
  // and finally the integration
  std::vector<int> localCubesIndexes = levels[thisLevel].getLocalCubesIndexes();
  const int N_local_cubes = localCubesIndexes.size();
  for (int i=0 ; i<N_local_cubes ; ++i) {
    int indexLocalCube = levels[thisLevel].cubesIndexesAfterReduction[localCubesIndexes[i]];
    levels[thisLevel].sphericalIntegration(ZI, levels[thisLevel].Sdown(indexLocalCube), levels[thisLevel].cubes[indexLocalCube], levels[thisLevel].thetas, levels[thisLevel].phis, w, mu_r, k, CFIE);
  }
}

void Octtree::computeFarField(blitz::Array<std::complex<float>, 2>& e_theta_far,
                              blitz::Array<std::complex<float>, 2>& e_phi_far,
                              const blitz::Array<float, 1>& r_phase_center,
                              const blitz::Array<float, 1>& octtreeXthetas_coarsest,
                              const blitz::Array<float, 1>& octtreeXphis_coarsest,
                              const blitz::Array<std::complex<float>, 1>& I_PQ,
                              const string octtree_data_path) /// coefficients of RWG functions
{
  blitz::Range all = blitz::Range::all();
  int ierror;
  if (this->getProcNumber()==0) cout << "\nFar field computation" << ": level ";

  const int N_levels = levels.size(), my_id = this->getProcNumber();
  int stopLevel;
  for (int l=0 ; l<N_levels ; ++l) {
    if (levels[l].DIRECTIONS_PARALLELIZATION==1) break;
    if (my_id==0) cout << levels[l].getLevel() << "..."; flush(cout);
    std::vector<int> localCubesIndexes(levels[l].getLocalCubesIndexes());
    const int N_local_cubes = localCubesIndexes.size();
    const int N_cubes = levels[l].getLevelSize();
    const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
    const int N_directions = (levels[l].DIRECTIONS_PARALLELIZATION!=1) ? N_theta*N_phi : levels[l].MPI_Scatterv_scounts(my_id);
    if (levels[l].Sdown.size()==0) levels[l].Sdown.resize(N_cubes);
    // we first compute the Sups of all the cubes at the given level

    if (l==0) { // we use computeSup only for the leaf cubes
      for (int i=0 ; i<N_local_cubes ; ++i) {
        int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
        if (levels[l].Sdown(indexLocalCube).size()==0) levels[l].Sdown(indexLocalCube).resize(2, N_theta*N_phi);
        levels[l].computeSup(levels[l].Sdown(indexLocalCube), k, I_PQ, levels[l].cubes[indexLocalCube], levels[l].thetas, levels[l].phis);
      }
    }
    else if ( l>0 ) { // the Sups are obtained by interpolation and shifting from the sonsCubes Sups
      blitz::Array<std::complex<float>, 2> S_tmp(2, N_theta*N_phi), S_tmp2(2, N_theta*N_phi);
      for (int i=0 ; i<N_local_cubes ; ++i) {
        int indexLocalCube = levels[l].cubesIndexesAfterReduction[localCubesIndexes[i]];
        if (levels[l].Sdown(indexLocalCube).size()==0) levels[l].Sdown(indexLocalCube).resize(2, N_directions);
        S_tmp2 = 0.0;
        const int NSons = levels[l].cubes[indexLocalCube].sonsIndexes.size();
        for (int j=0 ; j<NSons ; ++j) {
          const int sonIndex = levels[l-1].cubesIndexesAfterReduction[levels[l].cubes[indexLocalCube].sonsIndexes[j]];
          // interpolation
          interpolate2Dlfi(S_tmp(0, all), levels[l-1].Sdown(sonIndex)(0, all), levels[l-1].lfi2D);
          interpolate2Dlfi(S_tmp(1, all), levels[l-1].Sdown(sonIndex)(1, all), levels[l-1].lfi2D);
          // shifting
          const float * rc_1(levels[l].cubes[indexLocalCube].rCenter);
          const float * rc_2(levels[l-1].cubes[sonIndex].rCenter);
          const double DRcenters[3] = {rc_1[0]-rc_2[0], rc_1[1]-rc_2[1], rc_1[2]-rc_2[2]};
          shiftExp( S_tmp, levels[l].getShiftingArray(DRcenters[0], DRcenters[1], DRcenters[2]) );
          S_tmp2 += S_tmp;
        }
        levels[l].Sdown(indexLocalCube) = S_tmp2;
      }
    }
    stopLevel = l;
  }

  // we now define an interpolator from the last level we went through to the coarsest level.
  const int N_thetaCoarseLevel(octtreeXthetas_coarsest.size()), N_phiCoarseLevel(octtreeXphis_coarsest.size());
  // theta data
  float A_theta, B_theta;
  readFloatFromASCIIFile(octtree_data_path + "A_theta.txt", A_theta);
  readFloatFromASCIIFile(octtree_data_path + "B_theta.txt", B_theta);
  int INCLUDED_THETA_BOUNDARIES, PERIODIC_Theta, CYCLIC_Theta, NOrderInterpTheta;
  readIntFromASCIIFile(octtree_data_path + "INCLUDED_THETA_BOUNDARIES.txt", INCLUDED_THETA_BOUNDARIES);
  readIntFromASCIIFile(octtree_data_path + "PERIODIC_Theta.txt", PERIODIC_Theta);
  readIntFromASCIIFile(octtree_data_path + "CYCLIC_Theta.txt", CYCLIC_Theta);
  readIntFromASCIIFile(octtree_data_path + "NOrderInterpTheta.txt", NOrderInterpTheta);
  // phi data
  float A_phi, B_phi;
  readFloatFromASCIIFile(octtree_data_path + "A_phi.txt", A_phi);
  readFloatFromASCIIFile(octtree_data_path + "B_phi.txt", B_phi);
  int INCLUDED_PHI_BOUNDARIES, PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi;
  readIntFromASCIIFile(octtree_data_path + "INCLUDED_PHI_BOUNDARIES.txt", INCLUDED_PHI_BOUNDARIES);
  readIntFromASCIIFile(octtree_data_path + "PERIODIC_Phi.txt", PERIODIC_Phi);
  readIntFromASCIIFile(octtree_data_path + "CYCLIC_Phi.txt", CYCLIC_Phi);
  readIntFromASCIIFile(octtree_data_path + "NOrderInterpPhi.txt", NOrderInterpPhi);

  LagrangeFastInterpolator2D FarFieldInterpolator(octtreeXthetas_coarsest, levels[stopLevel].thetas, A_theta, B_theta, INCLUDED_THETA_BOUNDARIES, NOrderInterpTheta, PERIODIC_Theta, CYCLIC_Theta, octtreeXphis_coarsest, levels[stopLevel].phis, A_phi, B_phi, INCLUDED_PHI_BOUNDARIES, NOrderInterpPhi, PERIODIC_Phi, CYCLIC_Phi);

  // actual radiation pattern computation
  blitz::Array<std::complex<float>, 2> SupLastLevel(2, octtreeXthetas_coarsest.size() * octtreeXphis_coarsest.size()), SupLastLevelTmp(2, octtreeXthetas_coarsest.size() * octtreeXphis_coarsest.size());
  SupLastLevel = 0.0;
  std::vector<int> localCubesIndexes(levels[stopLevel].getLocalCubesIndexes());
  const int N_local_cubes = localCubesIndexes.size();
  for (int i=0 ; i<N_local_cubes ; ++i) {
    SupLastLevelTmp = 0.0;
    int indexLocalCube = levels[stopLevel].cubesIndexesAfterReduction[localCubesIndexes[i]];
    interpolate2Dlfi(SupLastLevelTmp(0, all), levels[stopLevel].Sdown(indexLocalCube)(0, all), FarFieldInterpolator);
    interpolate2Dlfi(SupLastLevelTmp(1, all), levels[stopLevel].Sdown(indexLocalCube)(1, all), FarFieldInterpolator);
    // shifting
    double DRcenters[3];
    for (int j=0 ; j<3 ; j++) DRcenters[j] = (r_phase_center(j) - levels[stopLevel].cubes[indexLocalCube].rCenter[j]);
    blitz::Array<std::complex<float>, 1> shiftingArray(octtreeXthetas_coarsest.size() * octtreeXphis_coarsest.size());
    for (int m=0 ; m<N_thetaCoarseLevel ; ++m) {
      const float sinTheta = sin(octtreeXthetas_coarsest(m));
      const float cosTheta = cos(octtreeXthetas_coarsest(m));
      for (int n=0 ; n<N_phiCoarseLevel ; ++n) {
        const float k_hat[3] = {sinTheta*cos(octtreeXphis_coarsest(n)), sinTheta*sin(octtreeXphis_coarsest(n)), cosTheta};
        shiftingArray(m + n*N_thetaCoarseLevel) = static_cast<std::complex<float> > (exp(-I*this->k * (k_hat[0]*DRcenters[0] + k_hat[1]*DRcenters[1] + k_hat[2]*DRcenters[2])));
      }
    }
    shiftExp( SupLastLevelTmp, shiftingArray );
    SupLastLevel += (static_cast<std::complex<float> >(-I*mu_0) * w/static_cast<float>(4.0*M_PI) * mu_r) * SupLastLevelTmp;
  }
  // we now gather all
  ierror = MPI_Barrier(MPI_COMM_WORLD);
  SupLastLevelTmp = 0.0;
  ierror = MPI_Allreduce(SupLastLevel.data(), SupLastLevelTmp.data(), SupLastLevelTmp.size(), MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  SupLastLevel = SupLastLevelTmp;
  e_theta_far.resize(N_thetaCoarseLevel, N_phiCoarseLevel);
  e_phi_far.resize(N_thetaCoarseLevel, N_phiCoarseLevel);
  // e_theta
  for (int m=0 ; m<N_thetaCoarseLevel ; ++m) {
    for (int n=0 ; n<N_phiCoarseLevel ; ++n) {
      e_theta_far(m, n) = SupLastLevel(0, m + n*N_thetaCoarseLevel);
    }
  }
  // e_phi
  for (int m=0 ; m<N_thetaCoarseLevel ; ++m) {
    for (int n=0 ; n<N_phiCoarseLevel ; ++n) {
      e_phi_far(m, n) = SupLastLevel(1, m + n*N_thetaCoarseLevel);
    }
  }
  if (this->getProcNumber()==0) {
    std::cout << "finished!" << std::endl;
    flush(std::cout);
  }
}

void Octtree::computeSourceFarField(blitz::Array<std::complex<float>, 2>& e_theta_far,
                                    blitz::Array<std::complex<float>, 2>& e_phi_far,
                                    const blitz::Array<float, 1>& r_phase_center,
                                    const blitz::Array<float, 1>& octtreeXthetas_coarsest,
                                    const blitz::Array<float, 1>& octtreeXphis_coarsest,
                                    const int J_DIPOLES_EXCITATION,
                                    const blitz::Array<std::complex<float>, 2>& J_dip,
                                    const blitz::Array<float, 2>& r_J_dip,
                                    const int M_DIPOLES_EXCITATION,
                                    const blitz::Array<std::complex<float>, 2>& M_dip,
                                    const blitz::Array<float, 2>& r_M_dip)
{
  const int my_id = this->getProcNumber();
  if (my_id==0) cout << "\nSource Far field computation...";
  const int N_thetas = octtreeXthetas_coarsest.size(), N_phis = octtreeXphis_coarsest.size();

  blitz::Array<std::complex<float>, 2> SupTmp(2, N_thetas*N_phis), Sup(2, N_thetas*N_phis);
  Sup = 0.0;
  if (J_DIPOLES_EXCITATION==1) {
    const int N_J_dipoles = J_dip.rows();
    for (int i=0; i<N_J_dipoles; i++) {
      const int IS_J_CURRENT = 1;
      const std::complex<float> J[3] = {J_dip(i, 0), J_dip(i, 1), J_dip(i, 2)};
      const float r_dip[3] = {r_J_dip(i, 0), r_J_dip(i, 1), r_J_dip(i, 2)};
      const float rCenter[3] = {r_phase_center(0), r_phase_center(1), r_phase_center(2)};
      computeDipoleSup(SupTmp, J, IS_J_CURRENT, r_dip, rCenter, octtreeXthetas_coarsest, octtreeXphis_coarsest);
      Sup += (static_cast<std::complex<float> >(-I*mu_0) * w/static_cast<float>(4.0*M_PI) * mu_r) * SupTmp;
    }
  }
  if (M_DIPOLES_EXCITATION==1) {
    const int N_M_dipoles = M_dip.rows();
    for (int i=0; i<N_M_dipoles; i++) {
      const int IS_J_CURRENT = 0;
      const std::complex<float> M[3] = {M_dip(i, 0), M_dip(i, 1), M_dip(i, 2)};
      const float r_dip[3] = {r_M_dip(i, 0), r_M_dip(i, 1), r_M_dip(i, 2)};
      const float rCenter[3] = {r_phase_center(0), r_phase_center(1), r_phase_center(2)};
      computeDipoleSup(SupTmp, M, IS_J_CURRENT, r_dip, rCenter, octtreeXthetas_coarsest, octtreeXphis_coarsest);
      Sup += static_cast<std::complex<float> >((I*k)/(4.0*M_PI)) * SupTmp;
    }
  }
  e_theta_far.resize(N_thetas, N_phis);
  e_phi_far.resize(N_thetas, N_phis);
  for (int m=0 ; m<N_thetas ; ++m) {
    for (int n=0 ; n<N_phis ; ++n) {
      e_theta_far(m, n) = Sup(0, m + n*N_thetas);
      e_phi_far(m, n) = Sup(1, m + n*N_thetas);
    }
  }
  if (my_id==0) {
    std::cout << "finished!" << std::endl;
    flush(std::cout);
  }
}

void Octtree::computeDipoleSup(blitz::Array<std::complex<float>, 2> & Sup,
                               const std::complex<float> J_dipole[3],
                               const int IS_J_CURRENT,
                               const float r_dipole[3],
                               const float rCenter[3],
                               const blitz::Array<float, 1>& thetas,
                               const blitz::Array<float, 1>& phis)
{
  const int NThetas = thetas.size(), NPhis = phis.size();

  blitz::Array< float, 2> kHats(NThetas * NPhis, 3), thetaHats(NThetas * NPhis, 3), phiHats(NThetas * NPhis, 3);
  blitz::Array< std::complex<float>, 2> FC3Components(NThetas * NPhis, 3);
  std::vector<float> sin_thetas, cos_thetas;
  sin_thetas.resize(NThetas);
  cos_thetas.resize(NThetas);
  for (int p=0 ; p<NThetas ; ++p) sincosf(thetas(p), &sin_thetas[p], &cos_thetas[p]);
  // initialisation of arrays
  for (int q=0 ; q<NPhis ; ++q) {
    const float cos_phi = cos(phis(q)), sin_phi = sin(phis(q));
    for (int p=0 ; p<NThetas ; ++p) {
      int index = p + q*NThetas;
      const float sin_theta = sin_thetas[p], cos_theta = cos_thetas[p];
      kHats(index, 0) = sin_theta*cos_phi;
      kHats(index, 1) = sin_theta*sin_phi;
      kHats(index, 2) = cos_theta;
      thetaHats(index, 0) = cos_theta*cos_phi;
      thetaHats(index, 1) = cos_theta*sin_phi;
      thetaHats(index, 2) = -sin_theta;
      phiHats(index, 0) = -sin_phi;
      phiHats(index, 1) = cos_phi;
      phiHats(index, 2) = 0.0;
      for (int i=0 ; i<3 ; i++) FC3Components(index, i) = 0.0;
    }
  }
  // computation of FC3Components array
  const std::complex<float> I_k(static_cast<std::complex<float> >(I*this->k));
  const std::complex<float> fj[3] = {J_dipole[0], J_dipole[1], J_dipole[2]};
  const float expArg[3] = {r_dipole[0]-rCenter[0], r_dipole[1]-rCenter[1], r_dipole[2]-rCenter[2]};
  for (int q=0 ; q<NPhis ; q++) {
    const int index_1 = q*NThetas;
    for (int p=0 ; p<NThetas ; p++) {
      const int index = p + index_1;
      const std::complex<float> a(I_k * (expArg[0]*kHats(index, 0) + expArg[1]*kHats(index, 1) + expArg[2]*kHats(index, 2))); 
      float c, s, e;
      e = (a.real() == 0.0) ? 1.0 : exp(a.real());
      sincosf(a.imag(), &s, &c);
      const std::complex<float> EXP(e * c, e * s);
      if (IS_J_CURRENT==1) { // we have electric dipole
        FC3Components(index, 0) += fj[0] * EXP;
        FC3Components(index, 1) += fj[1] * EXP;
        FC3Components(index, 2) += fj[2] * EXP;
      }
      else { // we have magnetic dipole
        FC3Components(index, 0) += (kHats(index, 1)*fj[2] - kHats(index, 2)*fj[1]) * EXP;
        FC3Components(index, 1) += (kHats(index, 2)*fj[0] - kHats(index, 0)*fj[2]) * EXP;
        FC3Components(index, 2) += (kHats(index, 0)*fj[1] - kHats(index, 1)*fj[0]) * EXP;
      }
    }
  } // end q loop
  // transformation from cartesian to spherical coordinates and assignation to Sup
  for (int i=0 ; i<Sup.extent(1) ; ++i) {
    Sup(0, i) = thetaHats(i, 0)*FC3Components(i, 0) + thetaHats(i, 1)*FC3Components(i, 1) + thetaHats(i, 2)*FC3Components(i, 2);
    Sup(1, i) = phiHats(i, 0)*FC3Components(i, 0) + phiHats(i, 1)*FC3Components(i, 1) + phiHats(i, 2)*FC3Components(i, 2);
  }
}

