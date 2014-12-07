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

#include "level.h"
#include "readWriteBlitzArrayFromFile.h"
#include "alpha_computation.h"


/****************************************************************************/
/********************************* Level ************************************/
/****************************************************************************/
Level::Level(){};
Level::Level(const int l,
             const double leaf_side_length,
             const double big_cube_lower_coord[3],
             const blitz::Array<double, 2>& cubes_centroids)
{
  numberTimesCopied = 0;
  level = l;
  DIRECTIONS_PARALLELIZATION = 0;
  cubeSideLength = leaf_side_length;
  maxNumberCubes1D = static_cast<int>(pow(2.0, level));
  int N_cubes_level_L = cubes_centroids.extent(0);
  cubes.reserve(N_cubes_level_L);
  for (int j=0 ; j<N_cubes_level_L ; ++j) {
    const double r_c[3] = {cubes_centroids(j, 0), cubes_centroids(j, 1), cubes_centroids(j, 2)};
    addNode(Cube(level, leaf_side_length, big_cube_lower_coord, r_c));
    cubes[j].setOldIndex(j); // the original index, from the python mesh
  }
}

Level::Level(const int l,
             const int N_expansion,
             const double leaf_side_length, 
             const double big_cube_lower_coord[3],
             const blitz::Array<double, 2>& cubes_centroids,
             const std::complex<double>& waveNumber,
             const blitz::Array<float, 1>& Xthetas,
             const blitz::Array<float, 1>& Wthetas,
             const int INCLUDED_THETA_BOUNDARIES,
             const int N_theta,
             const int PERIODIC_Theta,
             const int CYCLIC_Theta,
             const int NOrderInterpolatorTheta,
             const blitz::Array<float, 1>& Xphis,
             const blitz::Array<float, 1>& Wphis,
             const int INCLUDED_PHI_BOUNDARIES,
             const int N_phi,
             const int PERIODIC_Phi,
             const int CYCLIC_Phi,
             const int NOrderInterpolatorPhi,
             const blitz::Array<float, 1>& XthetasNextLevel,
             const blitz::Array<float, 1>& XphisNextLevel,
             const int VERBOSE) // leaf level constructor: no sons needed
{
  const int my_id = MPI::COMM_WORLD.Get_rank(), num_procs = MPI::COMM_WORLD.Get_size();
  numberTimesCopied = 0;
  k = waveNumber;
  level = l;
  leaf = true;
  ceiling = 0;
  DIRECTIONS_PARALLELIZATION = 0;
  N = N_expansion;
  cubeSideLength = leaf_side_length;
  maxNumberCubes1D = static_cast<int>(pow(2.0, level));
  if ( (my_id==0) && (VERBOSE==1) ) cout << "construction of the leaf (finest) level" << endl;
  int j, N_cubes_level_L = cubes_centroids.extent(0);
  cubes.resize(N_cubes_level_L);
  for (j=0 ; j<N_cubes_level_L ; ++j) {
    const double r_c[3] = {cubes_centroids(j, 0), cubes_centroids(j, 1), cubes_centroids(j, 2)};
    cubes[j] = Cube(level, leaf_side_length, big_cube_lower_coord, r_c);
    cubes[j].setOldIndex(j); // the original index, from the python mesh
  }
  if ( (my_id==0) && (VERBOSE==1) ) cout << "cubes.size() = " << cubes.size() << ", cubes.capacity() = " << cubes.capacity() << endl;
  sortCubesByParents();

  thetas.resize(N_theta);
  thetas = Xthetas;
  weightsThetas.resize(N_theta);
  weightsThetas = Wthetas;
  phis.resize(N_phi);
  phis = Xphis;
  weightsPhis.resize(N_phi);
  weightsPhis = Wphis;

  // Lagrange Fast Interpolator
  const float THETA_MIN = 0., THETA_MAX = M_PI;
  const float PHI_MIN = 0., PHI_MAX = 2.0*M_PI;
  lfi2D.setLfi2D(LagrangeFastInterpolator2D (XthetasNextLevel, Xthetas, THETA_MIN, THETA_MAX, INCLUDED_THETA_BOUNDARIES, NOrderInterpolatorTheta, PERIODIC_Theta, CYCLIC_Theta, XphisNextLevel, Xphis, PHI_MIN, PHI_MAX, INCLUDED_PHI_BOUNDARIES, NOrderInterpolatorPhi, PERIODIC_Phi, CYCLIC_Phi));

  if (my_id==0) cout << "The total number of directions at level " << getLevel() << " is N_theta*N_phi = " << N_theta << "*" << N_phi << " = " << N_theta * N_phi << endl << endl;
}

void Level::copyLevel(const Level & levelToCopy) // copy constructor
{
  const int my_id = MPI::COMM_WORLD.Get_rank();
  numberTimesCopied = levelToCopy.getNumberTimesCopied ();
  incrementNumberTimesCopied();
  leaf = levelToCopy.getLeaf();
  ceiling = levelToCopy.getCeiling();
  level = levelToCopy.getLevel();
  if (numberTimesCopied>15) cout << "    Warning: level " << level << " on Process " << my_id << " has been copied more than 15 times. numberTimesCopied = " << numberTimesCopied << endl;
  N = levelToCopy.getN();
  cubeSideLength = levelToCopy.getCubeSideLength();
  maxNumberCubes1D = levelToCopy.getMaxNumberCubes1D();
  NCubesX = levelToCopy.getNCubesX();
  NCubesY = levelToCopy.getNCubesY();
  NCubesZ = levelToCopy.getNCubesZ();
  offsetAlphaIndexX = levelToCopy.getOffsetAlphaIndexX();
  offsetAlphaIndexY = levelToCopy.getOffsetAlphaIndexY();
  offsetAlphaIndexZ = levelToCopy.getOffsetAlphaIndexZ();
  k = levelToCopy.getK();
  cubes.resize(levelToCopy.cubes.size());
  for (unsigned int i=0 ; i<levelToCopy.cubes.size() ; ++i) cubes[i] = levelToCopy.cubes[i];
  numbersToIndexes = levelToCopy.getNumbersToIndexes();
  localCubesIndexes = levelToCopy.getLocalCubesIndexes();
  listOfFcToBeReceived.resize(levelToCopy.getListOfFcToBeReceived().size());
  listOfFcToBeSent.resize(levelToCopy.getListOfFcToBeSent().size());
  listOfFcToBeReceived = levelToCopy.getListOfFcToBeReceived();
  listOfFcToBeSent = levelToCopy.getListOfFcToBeSent();
  cubesIndexesAfterReduction = levelToCopy.getCubesIndexesAfterReduction();
  const int N_theta = levelToCopy.getNThetas(), N_phi= levelToCopy.getNPhis();
  thetas.resize(N_theta);
  phis.resize(N_phi);
  thetas = levelToCopy.getThetas();
  phis = levelToCopy.getPhis();
  blitz::Array<int, 1> alphaTranslationsExtents(3);
  alphaTranslationsExtents = levelToCopy.getAlphaTranslationsExtents();
  const int Nx = alphaTranslationsExtents(0), Ny = alphaTranslationsExtents(1), Nz = alphaTranslationsExtents(2);

  alphaTranslations.resize(Nx, Ny, Nz);
  alphaTranslationsIndexesNonZeros.resize(Nx, Ny, Nz);
  for (int i=0 ; i<Nx ; ++i) {
    for (int j=0 ; j<Ny ; ++j) {
      for (int m=0 ; m<Nz ; ++m) {
        alphaTranslations(i,j,m).resize(levelToCopy.alphaTranslations(i,j,m).size());
        alphaTranslationsIndexesNonZeros(i,j,m).resize(levelToCopy.alphaTranslationsIndexesNonZeros(i,j,m).size());
        alphaTranslations(i,j,m) = levelToCopy.alphaTranslations(i,j,m);
        alphaTranslationsIndexesNonZeros(i,j,m) = levelToCopy.alphaTranslationsIndexesNonZeros(i,j,m);
      }
    }
  }
  alphaTranslationsIndexes.resize(alphaTranslationsIndexes.extent(0), alphaTranslationsIndexes.extent(1), alphaTranslationsIndexes.extent(2), alphaTranslationsIndexes.extent(3));
  alphaTranslationsIndexes = levelToCopy.alphaTranslationsIndexes;
  shiftingArrays.resize(levelToCopy.getShiftingArrays().extent(0), levelToCopy.getShiftingArrays().extent(1));
  shiftingArrays = levelToCopy.getShiftingArrays();
  weightsThetas.resize(levelToCopy.getWeightsThetas().size());
  weightsThetas = levelToCopy.getWeightsThetas();
  weightsPhis.resize(levelToCopy.getWeightsPhis().size());
  weightsPhis = levelToCopy.getWeightsPhis();
  lfi2D.setLfi2D(levelToCopy.getLfi2D());
  Sdown.resize(levelToCopy.Sdown.size());
  for (unsigned int i=0 ; i<Sdown.size() ; ++i) {
    Sdown(i).resize(levelToCopy.Sdown(i).extent(0), levelToCopy.Sdown(i).extent(1));
    Sdown(i) = levelToCopy.Sdown(i);
  }

  DIRECTIONS_PARALLELIZATION = levelToCopy.DIRECTIONS_PARALLELIZATION;
  MPI_Scatterv_scounts.resize(levelToCopy.MPI_Scatterv_scounts.size());
  MPI_Scatterv_displs.resize(levelToCopy.MPI_Scatterv_displs.size());
  MPI_Scatterv_scounts = levelToCopy.MPI_Scatterv_scounts;
  MPI_Scatterv_displs = levelToCopy.MPI_Scatterv_displs;
}

Level::Level(const Level& levelToCopy) // copy constructor
{
  copyLevel(levelToCopy);
}

Level& Level::operator=(const Level& levelToCopy) { // copy assignment
  copyLevel(levelToCopy);
  return *this;
}

Level::Level(const Level & sonLevel,
             const int N_expansion,
             const double big_cube_lower_coord[3],
             const blitz::Array<float, 1>& Xthetas,
             const blitz::Array<float, 1>& Wthetas,
             const int INCLUDED_THETA_BOUNDARIES,
             const int N_theta,
             const int PERIODIC_Theta,
             const int CYCLIC_Theta,
             const int NOrderInterpolatorTheta,
             const blitz::Array<float, 1>& Xphis,
             const blitz::Array<float, 1>& Wphis,
             const int INCLUDED_PHI_BOUNDARIES,
             const int N_phi,
             const int PERIODIC_Phi,
             const int CYCLIC_Phi,
             const int NOrderInterpolatorPhi,
             const blitz::Array<float, 1>& XthetasNextLevel,
             const blitz::Array<float, 1>& XphisNextLevel,
             const int VERBOSE)
{
  const int my_id = MPI::COMM_WORLD.Get_rank(), num_procs = MPI::COMM_WORLD.Get_size();
  numberTimesCopied = 0;
  N = N_expansion;
  k = sonLevel.getK();
  int j;
  level = sonLevel.getLevel()-1;
  leaf = false; // because this level is created from a lower (finer) level...
  ceiling = 0;
  DIRECTIONS_PARALLELIZATION = 0;
  cubeSideLength = 2.0*sonLevel.getCubeSideLength();
  maxNumberCubes1D = sonLevel.getMaxNumberCubes1D()/2;
  if ( (my_id==0) && (VERBOSE==1) ) cout << "construction of level " << level << endl;
  flush(cout);
  addNode(Cube(sonLevel.getCube(0), level, big_cube_lower_coord, cubeSideLength)); // initialization
  for (j=1 ; j<sonLevel.getLevelSize() ; ++j) // we walk through the sons list
  {
    if (sonLevel.getCube(j).getFatherNumber() == cubes.back().getNumber()) cubes.back().addSon(sonLevel.getCube(j));
    else this->addNode(Cube(sonLevel.getCube(j), level, big_cube_lower_coord, cubeSideLength));
  }
  vector<Cube>(cubes).swap(cubes); // swap trick for trimming exceeding capacity
  if ( (my_id==0) && (VERBOSE==1) ) cout << "cubes.size() = " << cubes.size() << ", cubes.capacity() = " << cubes.capacity() << endl;
  thetas.resize(N_theta);
  thetas = Xthetas;
  weightsThetas.resize(N_theta);
  weightsThetas = Wthetas;
  phis.resize(N_phi);
  phis = Xphis;
  weightsPhis.resize(N_phi);
  weightsPhis = Wphis;
  // Lagrange Fast Interpolator. The coarsest level does not need one...
  const float THETA_MIN = 0., THETA_MAX = M_PI;
  const float PHI_MIN = 0., PHI_MAX = 2.0*M_PI;
  if (level>2) lfi2D.setLfi2D(LagrangeFastInterpolator2D (XthetasNextLevel, Xthetas, THETA_MIN, THETA_MAX, INCLUDED_THETA_BOUNDARIES, NOrderInterpolatorTheta, PERIODIC_Theta, CYCLIC_Theta, XphisNextLevel, Xphis, PHI_MIN, PHI_MAX, INCLUDED_PHI_BOUNDARIES, NOrderInterpolatorPhi, PERIODIC_Phi, CYCLIC_Phi));
  // computation of the MPI_Scatterv_scounts / MPI_Scatterv_displs
  MPI_Scatterv_scounts.resize(num_procs);
  MPI_Scatterv_displs.resize(num_procs);
  const int N_directions = N_theta * N_phi;
  int displacement = 0;
  for (int i=0 ; i<num_procs ; ++i) {
    if (i<num_procs-1) MPI_Scatterv_scounts(i) = N_directions/num_procs;
    else MPI_Scatterv_scounts(i) = N_directions - displacement;
    MPI_Scatterv_displs(i) = displacement;
    displacement += MPI_Scatterv_scounts(i);
  }
  if ( (my_id==0) ) cout << "The total number of directions at level " << getLevel() << " is N_theta*N_phi = " << N_theta << "*" << N_phi << " = " << N_directions << endl << endl;
}

Level::~Level()
{
  cubes.clear();
  numbersToIndexes.clear();
  localCubesIndexes.clear();
  listOfFcToBeReceived.clear();
  listOfFcToBeSent.clear();
  cubesIndexesAfterReduction.clear();
  thetas.free();
  phis.free();
  weightsThetas.free();
  weightsPhis.free();
  alphaTranslations.free();
  alphaTranslationsIndexesNonZeros.free();
  alphaTranslationsIndexes.free();
  shiftingArrays.free();
  Sdown.free();
  MPI_Scatterv_scounts.free();
  MPI_Scatterv_displs.free();
}

void Level::NCubesXYZComputation(const int VERBOSE)
{
  int NxMax, NyMax, NzMax, NxMin, NyMin, NzMin, my_id = MPI::COMM_WORLD.Get_rank();
  NxMax = NyMax = NzMax = 0;
  NxMin = NyMin = NzMin = maxNumberCubes1D;
  for (unsigned int j=0 ; j<cubes.size() ; ++j) {
    const float * absoluteCartesianCoordTmp(cubes[j].absoluteCartesianCoord);
    if ( static_cast<int>(absoluteCartesianCoordTmp[0]) > NxMax ) NxMax = static_cast<int>(absoluteCartesianCoordTmp[0]);
    if ( static_cast<int>(absoluteCartesianCoordTmp[1]) > NyMax ) NyMax = static_cast<int>(absoluteCartesianCoordTmp[1]);
    if ( static_cast<int>(absoluteCartesianCoordTmp[2]) > NzMax ) NzMax = static_cast<int>(absoluteCartesianCoordTmp[2]);
    if ( static_cast<int>(absoluteCartesianCoordTmp[0]) < NxMin ) NxMin = static_cast<int>(absoluteCartesianCoordTmp[0]);
    if ( static_cast<int>(absoluteCartesianCoordTmp[1]) < NyMin ) NyMin = static_cast<int>(absoluteCartesianCoordTmp[1]);
    if ( static_cast<int>(absoluteCartesianCoordTmp[2]) < NzMin ) NzMin = static_cast<int>(absoluteCartesianCoordTmp[2]);
  }
  NCubesX = NxMax - NxMin + 1;
  NCubesY = NyMax - NyMin + 1;
  NCubesZ = NzMax - NzMin + 1;
  if ( (my_id==0) && (VERBOSE==1) ) cout << "Level " << this->level << " : NCubesX, NCubesY, NCubesZ = " <<  NCubesX << ", " << NCubesY << ", " << NCubesZ << endl;
  if (this->DIRECTIONS_PARALLELIZATION==1) {
    this->offsetAlphaIndexX = (this->getCeiling()) ? NCubesX-1 : min(NCubesX, 4)-1;
    this->offsetAlphaIndexY = (this->getCeiling()) ? NCubesY-1 : min(NCubesY, 4)-1;
    this->offsetAlphaIndexZ = (this->getCeiling()) ? NCubesZ-1 : min(NCubesZ, 4)-1;
  }
  else {
    offsetAlphaIndexX = 0;
    offsetAlphaIndexY = 0;
    offsetAlphaIndexZ = 0;
  }
  if ( (my_id==0) && (VERBOSE==1) ) cout << "Level " << this->level << " : offsetAlphaIndexX, offsetAlphaIndexY, offsetAlphaIndexZ = " << offsetAlphaIndexX << ", " << offsetAlphaIndexY << ", " << offsetAlphaIndexZ << endl;
}

void Level::alphaTranslationsComputation(const int VERBOSE,
                                         const float alphaTranslation_smoothing_factor,
                                         const float alphaTranslation_thresholdRelValueMax,
                                         const float alphaTranslation_RelativeCountAboveThreshold)
{
  Range all = Range::all();
  const int translationOrder = getN(), NThetas = getNThetas(), NPhis = getNPhis();
  const double lambda = static_cast<double>(2.0*M_PI/abs(getK()));
  int translationOrder_prime = static_cast<int>(ceil(translationOrder * alphaTranslation_smoothing_factor));
  const float cutting_coefficient = alphaTranslation_thresholdRelValueMax;
  int my_id = MPI::COMM_WORLD.Get_rank();
  int Nx = (this->getCeiling()) ? NCubesX : min(NCubesX, 4);
  int Ny = (this->getCeiling()) ? NCubesY : min(NCubesY, 4);
  int Nz = (this->getCeiling()) ? NCubesZ : min(NCubesZ, 4);
  if (this->DIRECTIONS_PARALLELIZATION==1) {
    // we cannot use the symmetries for the alpha matrices anymore
    Nx = 2*Nx-1;
    Ny = 2*Ny-1;
    Nz = 2*Nz-1;
  }

  blitz::Array<float, 2> thetasPhis_all_directions(NThetas * NPhis, 2), thetasPhis;
  blitz::Array<float, 1> weightsThetasPhis_all_directions(NThetas * NPhis), weightsThetasPhis;
  for (int i=0 ; i<NThetas ; ++i) {
    for (int j=0 ; j<NPhis ; ++j) {
      thetasPhis_all_directions(i + j*NThetas, 0) = thetas(i);
      thetasPhis_all_directions(i + j*NThetas, 1) = phis(j);
      weightsThetasPhis_all_directions(i + j*NThetas) = weightsThetas(i) * weightsPhis(j);
    }
  }
  int N_directions, N_total_directions = NThetas * NPhis;
  if ( this->DIRECTIONS_PARALLELIZATION==1 ) {
    N_directions = this->MPI_Scatterv_scounts(my_id);
    thetasPhis.resize(N_directions, 2);
    weightsThetasPhis.resize(N_directions);
    const int offset = this->MPI_Scatterv_displs(my_id);
    for (int i=0 ; i<N_directions ; ++i) {
      thetasPhis(i, all) = thetasPhis_all_directions(offset + i, all);
      weightsThetasPhis(i) = weightsThetasPhis_all_directions(offset + i);
    }
  }
  else {
    N_directions = N_total_directions;
    thetasPhis.resize(N_directions, 2);
    thetasPhis = thetasPhis_all_directions;
    weightsThetasPhis.resize(N_directions);
    weightsThetasPhis = weightsThetasPhis_all_directions;
  }
  thetasPhis_all_directions.free();
  weightsThetasPhis_all_directions.free();

  this->alphaTranslations.resize(Nx, Ny, Nz);
  this->alphaTranslationsIndexesNonZeros.resize(Nx, Ny, Nz);
  blitz::Array<std::complex<float>, 1> alpha(N_directions);
  //blitz::Array<std::complex<float>, 1> alpha_all_directions(N_total_directions);
  blitz::Array<int, 1> isAlphaNonZero(N_total_directions);
  if ( (my_id==0) && (VERBOSE==1) ) {
    cout << "\n    Process " << my_id << ", Level " << this->level << " alpha translations computation" << endl;
    cout << "    alpha.shape() = " << Nx << ", " << Ny << ", " << Nz << ", " << N_directions << endl;
    flush(cout);
  }

  double r_mn[3];
  for (int x = 0 ; x<Nx ; ++x) {
    if ( (my_id==0) && (VERBOSE==1) ) cout << "\r    " << (x+1)*100/Nx << " % computed";
    flush(cout);
    for (int y = 0 ; y<Ny ; ++y) {
      for (int z = 0 ; z<Nz ; ++z) {
        r_mn[0] = static_cast<double> (x-this->offsetAlphaIndexX);
        r_mn[1] = static_cast<double> (y-this->offsetAlphaIndexY);
        r_mn[2] = static_cast<double> (z-this->offsetAlphaIndexZ);
        if ( (abs(r_mn[0]) > 1.0) || (abs(r_mn[1]) > 1.0) || (abs(r_mn[2]) > 1.0) ) { /// if cartesian distance is sufficient
          r_mn[0] *= this->cubeSideLength;
          r_mn[1] *= this->cubeSideLength;
          r_mn[2] *= this->cubeSideLength;
          IT_theta_IT_phi_alpha_C2 (alpha, r_mn, getK(), translationOrder, translationOrder_prime, thetasPhis);
          alpha *= weightsThetasPhis;
          // we then seek the max of alpha_all_directions
          double max_abs_alpha_local = max(abs(alpha));
          double max_abs_alpha = max_abs_alpha_local;
          if ( this->DIRECTIONS_PARALLELIZATION==1 ) MPI_Allreduce(&max_abs_alpha_local, &max_abs_alpha, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          int countNonZeroLocal = 0, countNonZero;
          for (unsigned int kkk=0 ; kkk<alpha.size() ; ++kkk) {
            if (abs(alpha(kkk)) >= cutting_coefficient * max_abs_alpha) {
              countNonZeroLocal++;
              isAlphaNonZero(kkk) = 1;
            }
            else {
              alpha(kkk) = 0.0;
              isAlphaNonZero(kkk) = 0;
            }
          }
          /* To be erased */
          countNonZero = countNonZeroLocal;
          if ( this->DIRECTIONS_PARALLELIZATION==1 ) MPI_Allreduce(&countNonZeroLocal, &countNonZero, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
          if (countNonZero * 1.0/N_total_directions < alphaTranslation_RelativeCountAboveThreshold) {
            if ((my_id==0) && (this->cubeSideLength>= 16.0*lambda) && (VERBOSE==1) ) cout << "level " << this->getLevel() << ", cube side length = " << this->cubeSideLength/lambda << " lambdas, rel. count for r_mn = [" << r_mn[0]/this->cubeSideLength << ", " << r_mn[1]/this->cubeSideLength << ", " << r_mn[2]/this->cubeSideLength << "] is = " << countNonZero * 100.0/N_total_directions << endl;
            flush(cout);
            this->alphaTranslations(x, y, z).resize(countNonZeroLocal);
            this->alphaTranslationsIndexesNonZeros(x, y, z).resize(countNonZeroLocal);
            int index = 0;
            for (unsigned int kkk=0 ; kkk<alpha.size() ; ++kkk) {
              if (isAlphaNonZero(kkk) == 1) {
                this->alphaTranslations(x, y, z)(index) = alpha(kkk);
                this->alphaTranslationsIndexesNonZeros(x, y, z)(index) = kkk;
                index++;
              }
            }
            if (index != countNonZeroLocal) {
              cout << "error in computing alphaTranslations. index = " << index << ", countNonZero = " << countNonZero << ". Exiting..." << endl;
              exit(1);
            }
          }
          else {
            this->alphaTranslations(x, y, z).resize(alpha.size());
            this->alphaTranslations(x, y, z) = alpha;
            this->alphaTranslationsIndexesNonZeros(x, y, z).resize(0);
          }

        }
      }
    }
  }
  if ( (my_id==0) && (VERBOSE==1) ) cout << endl;
  // now we compute the alphaTranslationsIndexes. This is necessary due to the
  // fact that we use the symmetries in alpha computation, which results
  // in less computation time and less memory consumption for exactly the
  // same precision.
  if (this->DIRECTIONS_PARALLELIZATION!=1) {
    N_directions = NThetas * NPhis;
    alphaTranslationsIndexes.resize(2, 2, 2, N_directions);
    // if r_p < 0, then corresponding index in alphaTranslationsIndexes is 0.
    // for example, only r_y < 0: coordinates in alphaTranslationsIndexes are: (1, 0, 1)
    blitz::Array<int, 1> oldAlphaIndex(N_directions), newAlphaIndexZ(N_directions), newAlphaIndexY(N_directions), newAlphaIndexX(N_directions);
    for (int i=0; i<NThetas; ++i) {
      for (int j=0; j<NPhis; ++j) {
        const int index = i + j*NThetas;
        oldAlphaIndex(index) = index;
      }
    }
    for (int m=0; m<2; ++m) {
      for (int n=0; n<2; ++n) {
        for (int p=0; p<2; ++p) {
          alphaTranslationIndexConstructionX(newAlphaIndexX, oldAlphaIndex, m, NThetas, NPhis);
          alphaTranslationIndexConstructionY(newAlphaIndexY, newAlphaIndexX, n, NThetas, NPhis);
          alphaTranslationIndexConstructionZ(newAlphaIndexZ, newAlphaIndexY, p, NThetas, NPhis);
          alphaTranslationsIndexes(m, n, p, all) = newAlphaIndexZ;
        }
      }
    }
  }
}

double Level::getAlphaTranslationsSizeMB(void) const
{
  const int Nx = this->alphaTranslations.extent(0), Ny = this->alphaTranslations.extent(1), Nz = this->alphaTranslations.extent(2);
  int N_alpha_elements = 0;
  for (int x = 0 ; x<Nx ; ++x) {
    for (int y = 0 ; y<Ny ; ++y) {
      for (int z = 0 ; z<Nz ; ++z) {
        N_alpha_elements += this->alphaTranslations(x, y, z).size();
      }
    }
  }
  const double alphaTranslationsSizeMB = N_alpha_elements * 2.0 * 4.0 / (1024.0 * 1024.0);
  return alphaTranslationsSizeMB;
}

void Level::alphaTranslationIndexConstructionZ(blitz::Array<int, 1>& newAlphaIndex,
                                               const blitz::Array<int, 1>& oldAlphaIndex,
                                               const int alphaCartesianCoordZ,
                                               const int N_theta,
                                               const int N_phi)
/**
 * Symmetry with respect to the \f$ xOy \f$ plane of symmetry. We the have that
 * \f[
 * \alpha \left( \theta > \pi/2, \phi \right) = \alpha \left( \pi - \theta, \phi \right)
 * \f]
 * which, in indexes terms, can be translated as
 *
 */
{
  if (alphaCartesianCoordZ == 0) {
    for (int i=0; i<N_theta; ++i) {
      for (int j=0; j<N_phi; ++j) {
        const int index = i + j*N_theta;
        newAlphaIndex(index) = oldAlphaIndex(N_theta - i - 1 + j*N_theta);
      }
    }
  }
  else newAlphaIndex = oldAlphaIndex;
}
void Level::alphaTranslationIndexConstructionY(blitz::Array<int, 1>& newAlphaIndex,
                                               const blitz::Array<int, 1>& oldAlphaIndex,
                                               const int alphaCartesianCoordY,
                                               const int N_theta,
                                               const int N_phi)
/**
 * Symmetry with respect to the \f$ xOz \f$ plane of symmetry. We the have that
 * \f[
 * \alpha \left( \theta, \phi > \pi \right) = \alpha \left( \theta, 2 \pi - \phi \right)
 * \f]
 * which, in indexes terms, can be translated as
 *
 */
{
  if (alphaCartesianCoordY == 0) {
    for (int i=0; i<N_theta; ++i) {
      for (int j=0; j<N_phi; ++j) {
        const int index = i + j*N_theta;
        newAlphaIndex(index) = oldAlphaIndex(i + (N_phi-1-j) * N_theta);
      }
    }
  }
  else newAlphaIndex = oldAlphaIndex;
}
void Level::alphaTranslationIndexConstructionX(blitz::Array<int, 1>& newAlphaIndex,
                                               const blitz::Array<int, 1>& oldAlphaIndex,
                                               const int alphaCartesianCoordX,
                                               const int N_theta,
                                               const int N_phi)
/**
 * Symmetry with respect to the \f$ yOz \f$ plane of symmetry. We the have that
 * \f[
 * \alpha \left( \theta, \pi/2 < \phi < 3/2 \pi \right) = \alpha \left( \theta, \pi - \phi \right)
 * \f]
 * which, in indexes terms, can be translated as
 *
 */
{
  if (alphaCartesianCoordX == 0) {
    for (int i=0; i<N_theta; ++i) {
      for (int j=0; j<N_phi; ++j) {
        const int index = i + j*N_theta;
        if (j<N_phi/2) newAlphaIndex(index) = oldAlphaIndex(i + (N_phi/2-1-j) * N_theta);
        else newAlphaIndex(index) = oldAlphaIndex(i + (N_phi-1-(j-N_phi/2)) * N_theta);
      }
    }
  }
  else newAlphaIndex = oldAlphaIndex;
}

void Level::shiftingArraysComputation(void)
/**
 * This is the function that computes the shifting coefficients from each corner of the cube
 * towards its center. Therefore, there are as many shifting coefficients as there are elements
 * at the given level times 8, since the cube has 8 corners. 
 *
 * If one wants to go from the center towards a corner, as it happens in down-shifting, 
 * the coefficients must be conjugated.
 *
 * The "shiftingArrays" is in fact a 2D array, where the first column corresponds to the (0,0,0) corner
 * and the last column is the (1,1,1) corner. The numbering is as follows:
 * (0,0,0) -> column 0
 * (0,0,1) -> column 1
 * (0,1,0) -> column 2
 * (0,1,1) -> column 3
 * (1,0,0) -> column 4
 * (1,0,1) -> column 5
 * (1,1,0) -> column 6
 * (1,1,1) -> column 7
 */
{
  const int N_theta = thetas.size(), N_phi = phis.size();
  shiftingArrays.resize(8, N_theta * N_phi);
  double sinTheta, cosTheta, shiftingDistance = this->cubeSideLength/4.0;
  double k_hat[3], shiftingVector[3];
  blitz::Array< blitz::Array<double, 1>, 1> shiftingVectors(8);
  for (int p=0 ; p<8 ; ++p) shiftingVectors(p).resize(3);
  shiftingVectors(0) = -1.0, -1.0, -1.0;
  shiftingVectors(1) = -1.0, -1.0,  1.0;
  shiftingVectors(2) = -1.0,  1.0, -1.0;
  shiftingVectors(3) = -1.0,  1.0,  1.0;
  shiftingVectors(4) =  1.0, -1.0, -1.0;
  shiftingVectors(5) =  1.0, -1.0,  1.0;
  shiftingVectors(6) =  1.0,  1.0, -1.0;
  shiftingVectors(7) =  1.0,  1.0,  1.0;
  for (int p=0 ; p<8 ; ++p) {
    for (int i=0; i<3; i++) shiftingVector[i] = shiftingVectors(p)(i) * shiftingDistance;
    for (int m=0 ; m<N_theta ; ++m) {
      sinTheta = sin(thetas(m));
      cosTheta = cos(thetas(m));
      for (int n=0 ; n<N_phi ; ++n) {
        k_hat[0] = sinTheta*cos(phis(n));
        k_hat[1] = sinTheta*sin(phis(n));
        k_hat[2] = cosTheta;
        const double temp(k_hat[0]*shiftingVector[0] + k_hat[1]*shiftingVector[1] + k_hat[2]*shiftingVector[2]);
        shiftingArrays(p, m + n*N_theta) = static_cast<std::complex<float> > (exp(-I*k*temp));
      }
    }
  }
}

blitz::Array<int, 1> Level::getAlphaTranslationsExtents(void) const
{
  blitz::Array<int, 1> dimensions(3);
  dimensions = alphaTranslations.extent(0), alphaTranslations.extent(1), alphaTranslations.extent(2);
  return dimensions;
}

void Level::sortCubesByParents(void)
{
  sort(cubes.begin(), cubes.end());
  for (int j=0 ; j<getLevelSize() ; ++j) cubes[j].setIndex(j); /// we write the new values of the indexes
}

void Level::printCubesFathersNumbers(void) 
{
  int j, N = getLevelSize(), l = getLevel();
  for (j=0 ; j<N ; ++j) cout << "Level " << l << " : father number of cube " << j << " = " << cubes[j].getFatherNumber() << endl;
}

void Level::printCubesSonsIndexes(void) 
{
  int j, N = getLevelSize(), l = getLevel();
  vector<int> sonsIndexes;
  for (j=0 ; j<N ; ++j) {
    cout << "Level " << l << " : sons indexes of cube " << j << " = ";
    sonsIndexes = getCube(j).getSonsIndexes();
    for (unsigned int i=0; i<sonsIndexes.size(); ++i) cout << sonsIndexes[i] << " ";
    cout << endl;
  }
}

void Level::printCubesRCenters(void) 
{
  int j, N = getLevelSize(), l = getLevel();
  for (j=0 ; j<N ; ++j) cout << "Level " << l << " : center of cube " << j << " = " << cubes[j].rCenter[0] << ", " << cubes[j].rCenter[1] << ", " << cubes[j].rCenter[2] << ". Sidelength = " << getCubeSideLength() << endl;
}

void Level::printCubesRWG_numbers(void) 
{
  int N = getLevelSize(), l = getLevel();
  for (int j=0 ; j<N ; ++j) {
    std::cout << "Level " << l << " : RWG of cube " << j << " = ";
    for (unsigned int i=0 ; i<cubes[j].RWG_numbers.size() ; ++i) cout << cubes[j].RWG_numbers[i] << ", ";
    cout << endl;
  }
}

void Level::printThetaPhi(void) 
{
  cout << "thetas = " << thetas << endl;
  cout << "weightsThetas = " << weightsThetas << endl;
  cout << "phis = " << phis << endl;
  cout << "weightsPhis = " << weightsPhis << endl;
}

void Level::updateFatherIndexes(const Level& fatherLevel)
{
  vector<int> sonsIndexes;
  for (int i=0 ; i<fatherLevel.getLevelSize() ; ++i) {
    sonsIndexes = fatherLevel.getCube(i).getSonsIndexes();
    for (unsigned int j=0 ; j<sonsIndexes.size() ; ++j) {
      if (cubes[sonsIndexes[j]].getFatherNumber() != fatherLevel.getCube(i).getNumber()) {
        cout << "Level::updateFatherIndexes: father number in sons and father do not match!" << endl;
        exit(1);
      }
      cubes[sonsIndexes[j]].setFatherIndex(i);
    }
  }
}

void Level::setNumbersToIndexes(void)
{
  const int N_cubes = getLevelSize();
  int i, n;
  numbersToIndexes.reserve(N_cubes);
  for (i=0 ; i<N_cubes ; ++i) {
    n = getCube(i).getNumber();
    numbersToIndexes.push_back(Dictionary<int, int>(n, i));
  }
}

void Level::sortNumbersToIndexes(void)
{
  sort(numbersToIndexes.begin(), numbersToIndexes.end());
}

void Level::printNumbersToIndexes(void)
{
  const int N_cubes = getLevelSize();
  const int l = getLevel();
  for (int i=0 ; i<N_cubes ; ++i) {
    cout << "Level " << l << " : numbersToIndexes[" << i << "] = " << numbersToIndexes[i].getKey() << ", " << numbersToIndexes[i].getVal() << endl;
  }
}

void Level::printCubesNeighborsIndexes(void)
{
  int NCubes = getLevelSize(), l = getLevel();
  for (int j=0 ; j<NCubes ; ++j) {
    std::cout << "Level " << l << " : neighbors indexes of cube " << j << " = ";
    for (unsigned int i=0 ; i<getCubeNeighbors(j).size() ; ++i) cout << getCubeNeighbors(j)[i] << ", ";
    cout << endl;
  }
}

void Level::searchCubesNeighborsIndexes(void)
{
  setNumbersToIndexes();
  sortNumbersToIndexes();
  const int N_cubes = getLevelSize();
  for (int i=0 ; i<N_cubes ; ++i) {
    const float * absCartCoord(cubes[i].absoluteCartesianCoord);
    vector<int> neighborsIndexes;
    // we find the neighbors
    for (int x=-1 ; x<2 ; ++x) {
      for (int y=-1 ; y<2 ; ++y) {
        for (int z=-1 ; z<2 ; ++z) {
          int index = -1;
          const double CandidateAbsCartCoord[3] = {absCartCoord[0] + x, absCartCoord[1] + y, absCartCoord[2] + z};
          /// no component of (absoluteCartesianCoord(i) + p) -- where i=0,1,2 and p = x,y,z -- can be:
          /// (1) negative or (2) greater than MaxNumberCubes1D.
          int condition = 1;
          for (int j=0 ; j<3 ; ++j) condition *= ( (CandidateAbsCartCoord[j] >= 0) && (CandidateAbsCartCoord[j] < getMaxNumberCubes1D()) );
          /* we also do not want to consider the cube itself */
          condition *= !((x==0) && (y==0) && (z==0));
          if (condition>0) {
            int candidate_number = static_cast<int>( CandidateAbsCartCoord[0] * pow2(getMaxNumberCubes1D()) + CandidateAbsCartCoord[1] * getMaxNumberCubes1D() + CandidateAbsCartCoord[2] );
            index = getIndexOfNumber(candidate_number);
          }
          if (index>-1) neighborsIndexes.push_back(index);
        }
      }
    }
    // we now trim the excess capacity of neighborsIndexes
    cubes[i].neighborsIndexes.resize(neighborsIndexes.size());
    for (unsigned int j=0 ; j<neighborsIndexes.size() ; ++j) cubes[i].neighborsIndexes[j] = neighborsIndexes[j];
  }
}

int Level::getIndexOfNumber(const int number) const // "numbersToIndexes" must be ordered
                                                          // this is done in the calling function...
{
  int ind_inf, ind_sup, ind_mid, index;
  const int N = getNumbersToIndexesSize();
  if ( (number < getNumberToIndex(0)) || (number > getNumberToIndex(N-1)) ) index = -1;
  else {
    ind_inf = 0;
    ind_sup = N-1;
    while(ind_sup-ind_inf > 1) {
      ind_mid = (ind_sup+ind_inf)/2;
      if (number > getNumberToIndex(ind_mid)) ind_inf = ind_mid;
      else ind_sup = ind_mid;
    }
    if (number == getNumberToIndex(ind_inf)) index = getIndexToIndex(ind_inf);
    else if (number == getNumberToIndex(ind_sup)) index = getIndexToIndex(ind_sup);
    else index = -1;
  }
  return index;
}

void Level::computeLocalCubesIndexes(const int procNumber) {
  const int N_cubes = getLevelSize();
  localCubesIndexes.resize(0);
  if (this->DIRECTIONS_PARALLELIZATION==1) {
    for (int j=0 ; j<N_cubes ; ++j) localCubesIndexes.push_back(j);
  }
  else {
    for (int j=0 ; j<N_cubes ; ++j) {
      if (cubes[j].getProcNumber() == procNumber) localCubesIndexes.push_back(j);
    }
  }
  vector<int>(localCubesIndexes).swap(localCubesIndexes); // trick for trimming exceeding capacity
}

void Level::computeLevelReduction(void) {
  // this function aims at computing a new level, bookeeping only the 
  // (1) local cubes and (2) cubes whose Fc are to be received, so
  // that the original level can be copied and shrunk to a new level.
  const int NCubes = cubes.size();
  cubesIndexesAfterReduction.resize(NCubes);
  // initialization
  for (int i=0 ; i<NCubes ; ++i) cubesIndexesAfterReduction[i] = i; // this alone is no reduction!
  // now fill-in: we create a temporary list containing all the cubes we wanna keep
  vector<int> cubesToKeep;
  // first the really local cubes
  for (unsigned int i=0 ; i<localCubesIndexes.size() ; ++i) cubesToKeep.push_back(localCubesIndexes[i]);
  // then those whose Fc's are to be received
  for (unsigned int i=0 ; i<listOfFcToBeReceived.size() ; ++i) {
    vector<int> tmpList = listOfFcToBeReceived[i];
    for (unsigned int j=0 ; j<tmpList.size() ; ++j) {
      cubesToKeep.push_back(tmpList[j]);
    }
  }
  // now sorting
  sort(cubesToKeep.begin(), cubesToKeep.end());
  // now filling the interesting vectors!
  vector<Cube> newCubes;
  newCubes.resize(cubesToKeep.size());
  int newIndex = 0;
  for (unsigned int i=0 ; i<cubesToKeep.size() ; ++i) {
    int indexToKeep = cubesToKeep[i];
    cubesIndexesAfterReduction[indexToKeep] = newIndex;
    newIndex++;
    newCubes[i] = cubes[indexToKeep];
  }
  // now reduction!!
  cubes.swap(newCubes);
}

void Level::computeOldIndexesOfCubes(blitz::Array<int, 1>& oldIndexesOfCubes) {
  // this function computes the indexes of the cubes in the original mesh
  const int N_local_cubes = localCubesIndexes.size();
  oldIndexesOfCubes.resize(N_local_cubes);
  for (int i=0 ; i<N_local_cubes ; ++i) {
    int indexLocalCube = cubesIndexesAfterReduction[localCubesIndexes[i]];      
    oldIndexesOfCubes(i) = cubes[indexLocalCube].getOldIndex();
  }
}


void Level::computeGaussLocatedArguments(const blitz::Array<int, 1>& local_cubes_NRWG, 
                                         const blitz::Array<int, 1>& local_RWG_numbers, 
                                         const blitz::Array<int, 1>& local_RWG_Numbers_CFIE_OK, 
                                         const blitz::Array<int, 2>& local_RWGNumbers_signedTriangles,
                                         const blitz::Array<float, 2>& local_RWGNumbers_trianglesCoord,
                                         const int N_Gauss)
{
  if ( getLeaf() ) {
    const int N_local_cubes = localCubesIndexes.size();
    int startIndex_in_localArrays = 0;
    for (int i=0 ; i<N_local_cubes ; ++i) {
      const int indexLocalCube = cubesIndexesAfterReduction[localCubesIndexes[i]];      
      const int NRWG = local_cubes_NRWG(i);
      cubes[indexLocalCube].computeGaussLocatedArguments(local_RWG_numbers, local_RWG_Numbers_CFIE_OK, local_RWGNumbers_signedTriangles, local_RWGNumbers_trianglesCoord, startIndex_in_localArrays, NRWG, N_Gauss);
      startIndex_in_localArrays += NRWG;
    }
  } 
}

void Level::RWGs_renumbering(void)
{
  if ( getLeaf() ) {
    const int N_local_cubes = localCubesIndexes.size();
    int startIndex = 0;
    for (int i=0 ; i<N_local_cubes ; ++i) {
      const int indexLocalCube = cubesIndexesAfterReduction[localCubesIndexes[i]];      
      const int NRWG = cubes[indexLocalCube].RWG_numbers.size();
      for (int j=0; j<NRWG; j++) cubes[indexLocalCube].RWG_numbers[j] = startIndex + j;
      startIndex += NRWG;
    }
  } 
}

void Level::computeSup(blitz::Array<std::complex<float>, 2> & Sup,
                       const std::complex<double>& k,
                       const blitz::Array<std::complex<float>, 1>& I_PQ,
                       const Cube & cube,
                       const blitz::Array<float, 1>& thetas,
                       const blitz::Array<float, 1>& phis)
{
  const int NThetas = thetas.size(), NPhis = phis.size(), NGauss = cube.triangle_GaussCoord.extent(1)/3;

  // A. Francavilla (29-05-2013)
  double sum_weigths;
  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, NGauss);

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
  const std::complex<float> I_k(static_cast<std::complex<float> >(I*k));
  const float * rCenter = cube.rCenter;
  const int T = cube.Triangle_numberOfRWGs.size();
  int startIndex = 0, startIndex_r_opp = 0;
  for (int i=0; i<T; i++) {
    const int n_rwg = cube.Triangle_numberOfRWGs[i];
    for (int j=0; j<NGauss; j++) {
      const float r[3] = {cube.triangle_GaussCoord(i, j*3), cube.triangle_GaussCoord(i, j*3+1), cube.triangle_GaussCoord(i, j*3+2)};
      std::complex<float> fj[3] = {0.0, 0.0, 0.0};
      // loop on the RWGs for triangle i
      for (int rwg=0; rwg<n_rwg; rwg++) {
        const int RWG_index = cube.TriangleToRWGindex[startIndex + rwg];
        const float weight = cube.TriangleToRWGweight[startIndex + rwg] * weigths[j];
        const std::complex<float> i_pq = I_PQ(cube.RWG_numbers[RWG_index]) * weight;
        const int index = startIndex_r_opp + rwg*3;
        fj[0] += i_pq*(r[0]-cube.TriangleToRWG_ropp[index]);
        fj[1] += i_pq*(r[1]-cube.TriangleToRWG_ropp[index + 1]);
        fj[2] += i_pq*(r[2]-cube.TriangleToRWG_ropp[index + 2]);
      } // end loop RWGs
      const float expArg[3] = {r[0]-rCenter[0], r[1]-rCenter[1], r[2]-rCenter[2]};
      for (int q=0 ; q<NPhis/2 ; q++) {// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
        const int index_1 = q*NThetas, opp_index_1 = NThetas-1 + (q+NPhis/2) * NThetas;
        for (int p=0 ; p<NThetas ; p++) {
          const int index = p + index_1, opp_index = opp_index_1-p;
          const std::complex<float> a(I_k * (expArg[0]*kHats(index, 0) + expArg[1]*kHats(index, 1) + expArg[2]*kHats(index, 2))); 
          float c, s, e;
          e = (a.real() == 0.0) ? 1.0 : exp(a.real());
          sincosf(a.imag(), &s, &c);
          const std::complex<float> EXP(e * c, e * s);
          FC3Components(index, 0) += fj[0] * EXP;
          FC3Components(index, 1) += fj[1] * EXP;
          FC3Components(index, 2) += fj[2] * EXP;
          const std::complex<float> conjEXP((1.0/e) * c, -(1.0/e) * s);
          FC3Components(opp_index, 0) += fj[0] * conjEXP;
          FC3Components(opp_index, 1) += fj[1] * conjEXP;
          FC3Components(opp_index, 2) += fj[2] * conjEXP;
        }
      } // end q loop
    } // end j Gauss loop
    startIndex += n_rwg;
    startIndex_r_opp += n_rwg*3;
  }
  // transformation from cartesian to spherical coordinates and assignation to Sup
  for (int i=0 ; i<Sup.extent(1) ; ++i) {
    Sup(0, i) = thetaHats(i, 0)*FC3Components(i, 0) + thetaHats(i, 1)*FC3Components(i, 1) + thetaHats(i, 2)*FC3Components(i, 2);
    Sup(1, i) = phiHats(i, 0)*FC3Components(i, 0) + phiHats(i, 1)*FC3Components(i, 1) + phiHats(i, 2)*FC3Components(i, 2);
  }
}

void Level::sphericalIntegration(blitz::Array<std::complex<float>, 1>& ZI,
                                 const blitz::Array<std::complex<float>, 2>& Sdown,
                                 const Cube & cube,
                                 const blitz::Array<float, 1>& thetas,
                                 const blitz::Array<float, 1>& phis,
                                 const float w,
                                 const std::complex<float>& mu_r,
                                 const std::complex<double>& k,
                                 const blitz::Array<std::complex<float>, 1>& CFIE)
{
  const std::complex<float> ZZERO(0.0, 0.0);
  const bool nE_tmp = (CFIE(1)!=ZZERO); 
  const bool tH_tmp = (CFIE(2)!=ZZERO), nH_tmp = (CFIE(3)!=ZZERO);
  const std::complex<float> tJEFIE_factor(static_cast<std::complex<float> >(-I*mu_0)  * w * mu_r * CFIE(0)); // k²/(j*w*eps) factor
  const std::complex<float> nJEFIE_factor(static_cast<std::complex<float> >(-I*mu_0)  * w * mu_r * CFIE(1)); // k²/(j*w*eps) factor
  const std::complex<float> tJMFIE_factor(static_cast<std::complex<float> >(I*k) * CFIE(2)); // j*k factor
  const std::complex<float> nJMFIE_factor(static_cast<std::complex<float> >(I*k) * CFIE(3)); // j*k factor
  const int NThetas = thetas.size(), NPhis = phis.size(), NGauss = cube.triangle_GaussCoord.extent(1)/3;

  // A. Francavilla (29-05-2013)
  double sum_weigths;
  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, NGauss);

  blitz::Array< float, 2> kHats(NThetas * NPhis, 3);
  blitz::Array< std::complex<float>, 2> GC3Components(NThetas * NPhis, 3), GC3Exp(NThetas * NPhis, 3);
  std::vector<float> sin_thetas, cos_thetas;
  sin_thetas.resize(NThetas);
  cos_thetas.resize(NThetas);
  for (int p=0 ; p<NThetas ; ++p) sincosf(thetas(p), &sin_thetas[p], &cos_thetas[p]);
  // initialisation of arrays
  for (int q=0 ; q<NPhis ; ++q) {
    const float cos_phi = cos(phis(q)), sin_phi = sin(phis(q));
    const float phiHat[3] = {-sin_phi, cos_phi, 0.0};
    for (int p=0 ; p<NThetas ; ++p) {
      int index = p + q*NThetas;
      const float sin_theta = sin_thetas[p], cos_theta = cos_thetas[p];
      kHats(index, 0) = sin_theta*cos_phi;
      kHats(index, 1) = sin_theta*sin_phi;
      kHats(index, 2) = cos_theta;
      const float thetaHat[3] = {cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta};
      GC3Components(index, 0) = Sdown(0, index) * thetaHat[0] + Sdown(1, index) * phiHat[0];
      GC3Components(index, 1) = Sdown(0, index) * thetaHat[1] + Sdown(1, index) * phiHat[1];
      GC3Components(index, 2) = Sdown(0, index) * thetaHat[2] + Sdown(1, index) * phiHat[2];
    }
  }
  // computation of integration
  // defining local arrays used for faster computations
  const std::complex<float> minus_I_k(static_cast<std::complex<float> >(-I*k));
  const float * rCenter = cube.rCenter;
  const int T = cube.Triangle_numberOfRWGs.size();
  int startIndex = 0, startIndex_r_opp = 0;
  for (int i=0; i<T; i++) {
    const int n_rwg = cube.Triangle_numberOfRWGs[i];
    const float nHat[3] = {cube.triangle_nHat(i, 0), cube.triangle_nHat(i, 1), cube.triangle_nHat(i, 2)};
    for (int j=0; j<NGauss; j++) {
      const float r[3] = {cube.triangle_GaussCoord(i, j*3), cube.triangle_GaussCoord(i, j*3+1), cube.triangle_GaussCoord(i, j*3+2)};
      // computation of the shifting terms
      const float expArg[3] = {r[0]-rCenter[0], r[1]-rCenter[1], r[2]-rCenter[2]};
      std::complex<float> EJ[3] = {0.0, 0.0, 0.0};
      for (int q=0 ; q<NPhis/2 ; q++) {// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
        const int index_1 = q*NThetas, opp_index_1 = NThetas-1 + (q+NPhis/2) * NThetas;
        for (int p=0 ; p<NThetas ; p++) {
          const int index = p + index_1, opp_index = opp_index_1-p;
          const std::complex<float> a(minus_I_k * (expArg[0]*kHats(index, 0) + expArg[1]*kHats(index, 1) + expArg[2]*kHats(index, 2)));
          float c, s, e;
          e = (a.real() == 0.0) ? 1.0 : exp(a.real());
          sincosf(a.imag(), &s, &c);
          const std::complex<float> EXP(e * c, e * s);
          GC3Exp(index, 0) = GC3Components(index, 0)*EXP;
          GC3Exp(index, 1) = GC3Components(index, 1)*EXP;
          GC3Exp(index, 2) = GC3Components(index, 2)*EXP;
          const std::complex<float> conjEXP((1.0/e) * c, -(1.0/e) * s);
          GC3Exp(opp_index, 0) = GC3Components(opp_index, 0)*conjEXP;
          GC3Exp(opp_index, 1) = GC3Components(opp_index, 1)*conjEXP;
          GC3Exp(opp_index, 2) = GC3Components(opp_index, 2)*conjEXP;
          EJ[0] += (GC3Exp(index, 0) + GC3Exp(opp_index, 0));
          EJ[1] += (GC3Exp(index, 1) + GC3Exp(opp_index, 1));
          EJ[2] += (GC3Exp(index, 2) + GC3Exp(opp_index, 2));
        }
      }
      // now for the MFIEs
      std::complex<float> GC3_x_kHat[3] = {0.0, 0.0, 0.0};
      if (tH_tmp || nH_tmp) {
        // we first have to compute GC3_x_kHat
        for (int q=0 ; q<NPhis ; q++) {// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
          const int index_1 = q*NThetas;
          for (int p=0 ; p<NThetas ; p++) {
            const int index = p + index_1;
            GC3_x_kHat[0] += GC3Exp(index, 1) * kHats(index, 2) - GC3Exp(index, 2) * kHats(index, 1);
            GC3_x_kHat[1] += GC3Exp(index, 2) * kHats(index, 0) - GC3Exp(index, 0) * kHats(index, 2);
            GC3_x_kHat[2] += GC3Exp(index, 0) * kHats(index, 1) - GC3Exp(index, 1) * kHats(index, 0);
          }
        }
      } // end if for MFIE
      // loop on the RWGs for triangle i
      for (int rwg=0; rwg<n_rwg; rwg++) {
        // common EFIE and MFIE
        const int RWG_index = cube.TriangleToRWGindex[startIndex + rwg];
        const float weight = cube.TriangleToRWGweight[startIndex + rwg] * weigths[j];
        const int RWGNumber = cube.RWG_numbers[RWG_index];
        const int CFIE_OK = cube.RWG_numbers_CFIE_OK[RWG_index];
        // EFIE
        const int index = startIndex_r_opp + rwg*3;
        const float fj[3] = {(r[0]-cube.TriangleToRWG_ropp[index]), (r[1]-cube.TriangleToRWG_ropp[index + 1]), (r[2]-cube.TriangleToRWG_ropp[index + 2])};
        ZI(RWGNumber) += tJEFIE_factor * (EJ[0]*fj[0] + EJ[1]*fj[1] + EJ[2]*fj[2]) * weight;
        // we see if we need n x f_m
        const bool tH = tH_tmp * CFIE_OK;
        const bool nH = nH_tmp * CFIE_OK;
        float nHat_x_fj[3];
        if (nE_tmp || nH) {
          nHat_x_fj[0] = nHat[1]*fj[2]-nHat[2]*fj[1];
          nHat_x_fj[1] = nHat[2]*fj[0]-nHat[0]*fj[2];
          nHat_x_fj[2] = nHat[0]*fj[1]-nHat[1]*fj[0];
        }
        if (nE_tmp) ZI(RWGNumber) += nJEFIE_factor * (EJ[0]*nHat_x_fj[0] + EJ[1]*nHat_x_fj[1] + EJ[2]*nHat_x_fj[2]) * weight;
        // MFIE
        if (tH) ZI(RWGNumber) += tJMFIE_factor * (fj[0]*GC3_x_kHat[0] + fj[1]*GC3_x_kHat[1] + fj[2]*GC3_x_kHat[2]) * weight;
        if (nH) ZI(RWGNumber) += nJMFIE_factor * (nHat_x_fj[0]*GC3_x_kHat[0] + nHat_x_fj[1]*GC3_x_kHat[1] + nHat_x_fj[2]*GC3_x_kHat[2]) * weight;
      }// end loop on the RWGs
    }// end Gauss integration points loop
    startIndex += n_rwg;
    startIndex_r_opp += n_rwg*3;
  }// end triangles loop
}

