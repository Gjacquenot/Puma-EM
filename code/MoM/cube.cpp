#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace blitz;

#include "cube.h"
#include "readWriteBlitzArrayFromFile.h"

Cube::Cube(const bool is_leaf,                           // 1 if cube is leaf
           const int level,                              // the level
           const double sideLength,                      // length of cube side
           const blitz::TinyVector<double, 3>& bigCubeLowerCoord, // coordinates of level 0 cube
           const blitz::Array<double, 1>& r_c)                  // coordinates of center
{
  leaf = is_leaf;
  for (int i=0 ; i<3 ; ++i) rCenter(i) = r_c(i); // we must loop, since rCenter is a TinyVector

  // we compute the absolute cartesian coordinates and the cube number
  absoluteCartesianCoord = floor( (rCenter-bigCubeLowerCoord)/sideLength );
  double maxNumberCubes1D = pow(2.0, level);
  number = static_cast<int>( absoluteCartesianCoord(0) * pow2(maxNumberCubes1D) + absoluteCartesianCoord(1) * maxNumberCubes1D + absoluteCartesianCoord(2) );

  // we compute the number of the father
  blitz::TinyVector<double, 3> cartesianCoordInFathers = floor( (rCenter-bigCubeLowerCoord)/(2.0*sideLength) );
  double maxNumberCubes1D_next_level = maxNumberCubes1D/2.0;
  fatherNumber =  static_cast<int>( cartesianCoordInFathers(0) * pow2(maxNumberCubes1D_next_level) + cartesianCoordInFathers(1) * maxNumberCubes1D_next_level + cartesianCoordInFathers(2) );
}

Cube::Cube(const Cube& sonCube,
           const int level,
           const blitz::TinyVector<double, 3>& bigCubeLowerCoord,
           const double sideLength)
{
  leaf = 0; // since we construct from a son cube...
  number = sonCube.getFatherNumber();
  procNumber = sonCube.getProcNumber();
  sonsIndexes.push_back(sonCube.getIndex());
  blitz::TinyVector<double, 3> sonCartesianCoordInFathers = floor( (sonCube.getRCenter() - bigCubeLowerCoord) / sideLength );
  rCenter = bigCubeLowerCoord + sonCartesianCoordInFathers * sideLength + sideLength/2.0;
  // we compute the absolute cartesian coordinates
  absoluteCartesianCoord = floor( (rCenter-bigCubeLowerCoord)/sideLength );

  // we compute the number of the father of _this_ cube
  // (i.e. grandfather of sonCube)
  blitz::TinyVector<double, 3> cartesianCoordInFathers = floor( (rCenter-bigCubeLowerCoord)/(2.0*sideLength) );
  double maxNumberCubes1D_next_level = pow(2.0, level-1);
  fatherNumber = static_cast<int>(cartesianCoordInFathers(0) * pow2(maxNumberCubes1D_next_level) + cartesianCoordInFathers(1) * maxNumberCubes1D_next_level + cartesianCoordInFathers(2));
}

void Cube::computeGaussLocatedArguments(const blitz::Array<int, 1>& local_RWG_numbers,
                                        const blitz::Array<int, 1>& local_RWG_Numbers_CFIE_OK,
                                        const blitz::Array<float, 2>& local_RWGNumbers_trianglesCoord,
                                        const int startIndex_in_localArrays,
                                        const int NRWG,
                                        const int N_Gauss)
{
  RWG_numbers.resize(NRWG);
  RWG_numbers_CFIE_OK.resize(NRWG);
  for (int j=0 ; j<NRWG ; ++j) RWG_numbers[j] = local_RWG_numbers(startIndex_in_localArrays + j);
  for (int j=0 ; j<NRWG ; ++j) RWG_numbers_CFIE_OK[j] = local_RWG_Numbers_CFIE_OK(startIndex_in_localArrays + j);
  
  GaussLocatedWeightedRWG.resize(NRWG, 2*N_Gauss);
  GaussLocatedWeighted_nHat_X_RWG.resize(NRWG, 2*N_Gauss);
  GaussLocatedExpArg.resize(NRWG, 2*N_Gauss);
  double sum_weigths;
  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_Gauss);

  for (int j=0 ; j<NRWG ; ++j) {
    blitz::TinyVector<double, 3> r, r_p, r0, r1, r2, n_hat;
    for (int halfBasisCounter = 0 ; halfBasisCounter < 2 ; ++halfBasisCounter) {
      if (halfBasisCounter==0) {
        for (int i=0; i<3; i++) {
          r0(i) = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + j, i);
          r1(i) = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + j, i+3);
          r2(i) = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + j, i+6);
        }
      }
      else {
        for (int i=0; i<3; i++) {
          r0(i) = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + j, i+9);
          r1(i) = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + j, i+6);
          r2(i) = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + j, i+3);
        }
      }
      r_p = r0;
      const blitz::TinyVector<double, 3> r1_r0(r1 - r0), r2_r0(r2-r0);
      n_hat = cross(r1_r0, r2_r0);
      const double Area = sqrt(dot(n_hat, n_hat))/2.0;
      n_hat = n_hat/(2.0*Area);
      double l_p = sqrt(dot((r2-r1), (r2-r1)));
      double sign_edge_p = (halfBasisCounter==0) ? 1.0 : -1.0;
      for (int i=0 ; i<N_Gauss ; ++i) {
        r = r0 * xi[i] + r1 * eta[i] + r2 * (1-xi[i]-eta[i]);
        const blitz::TinyVector<double, 3> r_rp(r - r_p);
        GaussLocatedWeightedRWG(j, i + halfBasisCounter*N_Gauss) = (sign_edge_p * l_p/2.0/sum_weigths * weigths[i]) * (r_rp);
        GaussLocatedWeighted_nHat_X_RWG(j, i + halfBasisCounter*N_Gauss) = (sign_edge_p * l_p/2.0/sum_weigths * weigths[i]) * cross(n_hat, r_rp);
        GaussLocatedExpArg(j, i + halfBasisCounter*N_Gauss) = (r-rCenter);
      }
    }
  }
}


void Cube::copyCube(const Cube& cubeToCopy) // copy member function
{
  leaf = cubeToCopy.getLeaf();
  number = cubeToCopy.getNumber();
  index = cubeToCopy.getIndex();
  oldIndex = cubeToCopy.getOldIndex();
  procNumber = cubeToCopy.getProcNumber();
  fatherNumber = cubeToCopy.getFatherNumber();
  fatherProcNumber = cubeToCopy.getFatherProcNumber();
  fatherIndex = cubeToCopy.getFatherIndex();
  sonsIndexes = cubeToCopy.getSonsIndexes();
  sonsProcNumbers = cubeToCopy.getSonsProcNumbers();
  neighborsIndexes.resize(cubeToCopy.neighborsIndexes.size());
  neighborsIndexes = cubeToCopy.neighborsIndexes;
  localAlphaTransParticipantsIndexes.resize(cubeToCopy.localAlphaTransParticipantsIndexes.size());
  localAlphaTransParticipantsIndexes = cubeToCopy.localAlphaTransParticipantsIndexes;
  nonLocalAlphaTransParticipantsIndexes.resize(cubeToCopy.nonLocalAlphaTransParticipantsIndexes.size());
  nonLocalAlphaTransParticipantsIndexes = cubeToCopy.nonLocalAlphaTransParticipantsIndexes;
  rCenter = cubeToCopy.getRCenter();
  absoluteCartesianCoord = cubeToCopy.getAbsoluteCartesianCoord();
  RWG_numbers.resize(cubeToCopy.RWG_numbers.size());
  RWG_numbers = cubeToCopy.RWG_numbers;
  RWG_numbers_CFIE_OK.resize(cubeToCopy.RWG_numbers_CFIE_OK.size());
  RWG_numbers_CFIE_OK = cubeToCopy.RWG_numbers_CFIE_OK;
  GaussLocatedWeightedRWG.resize(cubeToCopy.GaussLocatedWeightedRWG.extent(0), cubeToCopy.GaussLocatedWeightedRWG.extent(1));
  GaussLocatedWeighted_nHat_X_RWG.resize(cubeToCopy.GaussLocatedWeighted_nHat_X_RWG.extent(0), cubeToCopy.GaussLocatedWeighted_nHat_X_RWG.extent(1));
  GaussLocatedExpArg.resize(cubeToCopy.GaussLocatedExpArg.extent(0), cubeToCopy.GaussLocatedExpArg.extent(1));
  GaussLocatedWeightedRWG = cubeToCopy.getGaussLocatedWeightedRWG();
  GaussLocatedWeighted_nHat_X_RWG = cubeToCopy.getGaussLocatedWeighted_nHat_X_RWG();
  GaussLocatedExpArg = cubeToCopy.getGaussLocatedExpArg();
}

Cube::Cube(const Cube& cubeToCopy) // copy constructor
{
  copyCube(cubeToCopy);
}

Cube& Cube::operator=(const Cube& cubeToCopy) { // copy assignment
  copyCube(cubeToCopy);
  return *this;
}

Cube::~Cube() {
  sonsIndexes.clear();
  sonsProcNumbers.clear();
  neighborsIndexes.clear();
  localAlphaTransParticipantsIndexes.clear();
  nonLocalAlphaTransParticipantsIndexes.clear();
  RWG_numbers.clear();
  RWG_numbers_CFIE_OK.clear();
  GaussLocatedWeightedRWG.free();
  GaussLocatedWeighted_nHat_X_RWG.free();
  GaussLocatedExpArg.free();
}

void Cube::addSon(const Cube& sonCube)
{
  if (number == sonCube.getFatherNumber()) sonsIndexes.push_back(sonCube.getIndex());
  else {
    cout << "ERROR: no Son added because (number != sonCube.getFatherNumber())" << endl;
    exit(1);
  }
}


bool Cube::operator== (const Cube & right) const {
  if ( this->getFatherNumber() == right.getFatherNumber() ) return 1;
  else return 0;
}

bool Cube::operator< (const Cube & right) const {
  if ( this->getFatherNumber() < right.getFatherNumber() ) return 1;
  else return 0;
}

