#include <mpi.h>

using namespace std;

#include "cube.h"
#include "readWriteBlitzArrayFromFile.h"

Cube::Cube(const int level,                              // the level
           const double sideLength,                      // length of cube side
           const double bigCubeLowerCoord[3], // coordinates of level 0 cube
           const double r_c[3])                  // coordinates of center
{
  for (int i=0 ; i<3 ; ++i) rCenter[i] = r_c[i]; // we must loop, since rCenter is an array

  // we compute the absolute cartesian coordinates and the cube number
  for (int i=0 ; i<3 ; ++i) absoluteCartesianCoord[i] = floor( (rCenter[i]-bigCubeLowerCoord[i])/sideLength );
  double maxNumberCubes1D = pow(2.0, level);
  number = static_cast<int>( absoluteCartesianCoord[0] * maxNumberCubes1D*maxNumberCubes1D + absoluteCartesianCoord[1] * maxNumberCubes1D + absoluteCartesianCoord[2] );

  // we compute the number of the father
  double cartesianCoordInFathers[3];
  for (int i=0; i<3; i++) cartesianCoordInFathers[i] = floor( (rCenter[i]-bigCubeLowerCoord[i])/(2.0*sideLength) );
  double maxNumberCubes1D_next_level = maxNumberCubes1D/2.0;
  fatherNumber =  static_cast<int>( cartesianCoordInFathers[0] * maxNumberCubes1D_next_level*maxNumberCubes1D_next_level + cartesianCoordInFathers[1] * maxNumberCubes1D_next_level + cartesianCoordInFathers[2] );
  // resizing the arrays to 0
  RWG_numbers.resize(0);
  RWG_numbers_CFIE_OK.resize(0);
  Triangle_numberOfRWGs.resize(0);
  TriangleToRWGindex.resize(0);
  TriangleToRWGweight.resize(0); 
  TriangleToRWG_ropp.resize(0);
  triangle_GaussCoord.resize(0, 0);
  triangle_nHat.resize(0, 0);
}

Cube::Cube(const Cube& sonCube,
           const int level,
           const double bigCubeLowerCoord[3],
           const double sideLength)
{
  number = sonCube.getFatherNumber();
  procNumber = sonCube.getProcNumber();
  sonsIndexes.push_back(sonCube.getIndex());
  double sonCartesianCoordInFathers[3];
  for (int i=0; i<3; i++) sonCartesianCoordInFathers[i] = floor( (sonCube.rCenter[i] - bigCubeLowerCoord[i]) / sideLength );
  for (int i=0; i<3; i++) rCenter[i] = bigCubeLowerCoord[i] + sonCartesianCoordInFathers[i] * sideLength + sideLength/2.0;
  // we compute the absolute cartesian coordinates
  for (int i=0 ; i<3 ; ++i) absoluteCartesianCoord[i] = floor( (rCenter[i]-bigCubeLowerCoord[i])/sideLength );

  // we compute the number of the father of _this_ cube
  // (i.e. grandfather of sonCube)
  double cartesianCoordInFathers[3];
  for (int i=0; i<3; i++) cartesianCoordInFathers[i] = floor( (rCenter[i]-bigCubeLowerCoord[i])/(2.0*sideLength) );
  double maxNumberCubes1D_next_level = pow(2.0, level-1);
  fatherNumber = static_cast<int>(cartesianCoordInFathers[0] * maxNumberCubes1D_next_level*maxNumberCubes1D_next_level + cartesianCoordInFathers[1] * maxNumberCubes1D_next_level + cartesianCoordInFathers[2]);
  // resizing the arrays to 0
  RWG_numbers.resize(0);
  RWG_numbers_CFIE_OK.resize(0);
  Triangle_numberOfRWGs.resize(0);
  TriangleToRWGindex.resize(0);
  TriangleToRWGweight.resize(0); 
  TriangleToRWG_ropp.resize(0);
  triangle_GaussCoord.resize(0, 0);
  triangle_nHat.resize(0, 0);
}

void Cube::computeGaussLocatedArguments(const blitz::Array<int, 1>& local_RWG_numbers,
                                        const blitz::Array<int, 1>& local_RWG_Numbers_CFIE_OK,
                                        const blitz::Array<int, 2>& local_RWGNumbers_signedTriangles,
                                        const blitz::Array<float, 2>& local_RWGNumbers_trianglesCoord,
                                        const int startIndex_in_localArrays,
                                        const int NRWG,
                                        const int N_Gauss)
{
  RWG_numbers.resize(NRWG);
  RWG_numbers_CFIE_OK.resize(NRWG);
  for (int j=0 ; j<NRWG ; ++j) RWG_numbers[j] = local_RWG_numbers(startIndex_in_localArrays + j);
  for (int j=0 ; j<NRWG ; ++j) RWG_numbers_CFIE_OK[j] = local_RWG_Numbers_CFIE_OK(startIndex_in_localArrays + j);

  double sum_weigths;
  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_Gauss);
  
  std::vector< Dictionary2<int, int, float> > DictTriangleToRWG;
  DictTriangleToRWG.reserve(NRWG*2);

  for (int j=0 ; j<NRWG ; j++) {
    const int triangle_number1 = abs(local_RWGNumbers_signedTriangles(startIndex_in_localArrays + j, 0));
    const int triangle_number2 = abs(local_RWGNumbers_signedTriangles(startIndex_in_localArrays + j, 1));
    DictTriangleToRWG.push_back(Dictionary2<int, int, float> (triangle_number1, j, 1.0));
    DictTriangleToRWG.push_back(Dictionary2<int, int, float> (triangle_number2, j, -1.0));
  }
  sort(DictTriangleToRWG.begin(), DictTriangleToRWG.end());
  // we count the number of triangles
  unsigned int T = 1;
  for (unsigned int j=1; j<DictTriangleToRWG.size(); j++) {
    if (DictTriangleToRWG[j].getKey() != DictTriangleToRWG[j-1].getKey()) T++;
  }
  // we now construct a matrix linking triangles and RWGs
  std::vector< std::vector<int> > TriangleToRWGindexTmp;
  std::vector< std::vector<float> > TriangleToRWGweightTmp, TriangleToRWG_roppTmp;
  TriangleToRWGindexTmp.reserve(T);
  TriangleToRWGweightTmp.reserve(T);
  TriangleToRWG_roppTmp.reserve(T);
  std::vector<int> index_tmp;
  std::vector<float> sign_tmp;
  // init
  index_tmp.push_back(DictTriangleToRWG[0].getVal1());
  sign_tmp.push_back(DictTriangleToRWG[0].getVal2());
  // loop
  for (unsigned int j=1; j<DictTriangleToRWG.size(); j++) {
    if (DictTriangleToRWG[j].getKey() != DictTriangleToRWG[j-1].getKey()) {
      std::vector<int>(index_tmp).swap(index_tmp);
      std::vector<float>(sign_tmp).swap(sign_tmp);
      TriangleToRWGindexTmp.push_back(index_tmp);
      TriangleToRWGweightTmp.push_back(sign_tmp);
      index_tmp.resize(0);
      sign_tmp.resize(0);
    }
    index_tmp.push_back(DictTriangleToRWG[j].getVal1());
    sign_tmp.push_back(DictTriangleToRWG[j].getVal2());
  }
  if (TriangleToRWGindex.size() < T) {
    if (index_tmp.size()!=0) {
      std::vector<int>(index_tmp).swap(index_tmp);
      std::vector<float>(sign_tmp).swap(sign_tmp);
      TriangleToRWGindexTmp.push_back(index_tmp);
      TriangleToRWGweightTmp.push_back(sign_tmp);
    }
    else {
      std::cout << "Error in the construction of cubes. Aborting" << std::endl;
      exit(1);
    }
  } 

  // construction of triangle_GaussCoord
  triangle_GaussCoord.resize(T, N_Gauss*3);
  triangle_nHat.resize(T, 3);
  for (unsigned int j=0; j<T; j++) {
    const int RWG_index = TriangleToRWGindexTmp[j][0];
    const float sign = TriangleToRWGweightTmp[j][0];
    double r[3], r0[3], r1[3], r2[3], n_hat[3], r1_r0[3], r2_r0[3];
    if (sign>0.0) {
      for (int i=0; i<3; i++) {
        r0[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i);
        r1[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+3);
        r2[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+6);
        r1_r0[i] = r1[i] - r0[i];
        r2_r0[i] = r2[i] - r0[i];
      }
    }
    else {
      for (int i=0; i<3; i++) {
        r0[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+9);
        r1[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+6);
        r2[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+3);
        r1_r0[i] = r1[i] - r0[i];
        r2_r0[i] = r2[i] - r0[i];
      }
    }
    // triangle normal computation
    n_hat[0] = r1_r0[1]*r2_r0[2] - r1_r0[2]*r2_r0[1];
    n_hat[1] = r1_r0[2]*r2_r0[0] - r1_r0[0]*r2_r0[2];
    n_hat[2] = r1_r0[0]*r2_r0[1] - r1_r0[1]*r2_r0[0];
    const double Area = 0.5*sqrt(n_hat[0]*n_hat[0] + n_hat[1]*n_hat[1] + n_hat[2]*n_hat[2]);
    for (int i=0; i<3; i++) triangle_nHat(int(j), i) = n_hat[i] * 1.0/(2.0*Area);
    // Gauss coord in triangles computation
    for (int i=0 ; i<N_Gauss ; ++i) {
      r[0] = r0[0] * xi[i] + r1[0] * eta[i] + r2[0] * (1.0-xi[i]-eta[i]);
      r[1] = r0[1] * xi[i] + r1[1] * eta[i] + r2[1] * (1.0-xi[i]-eta[i]);
      r[2] = r0[2] * xi[i] + r1[2] * eta[i] + r2[2] * (1.0-xi[i]-eta[i]);
      triangle_GaussCoord(int(j), i*3) = r[0];
      triangle_GaussCoord(int(j), i*3+1) = r[1];
      triangle_GaussCoord(int(j), i*3+2) = r[2];
    }
  }

  // computation of the weights and opposite vector for each RWG
  Triangle_numberOfRWGs.resize(T);
  int N_total_RWG = 0;
  for (unsigned int j=0; j<T; j++) {
    const int n_rwg = TriangleToRWGindexTmp[j].size();
    Triangle_numberOfRWGs[j] = n_rwg;
    N_total_RWG += n_rwg;
    std::vector<float> r_p;
    for (int p=0; p<n_rwg; p++) {
      const int RWG_index = TriangleToRWGindexTmp[j][p];
      const float sign = TriangleToRWGweightTmp[j][p];
      double r0[3], r1[3], r2[3], r2_r1[3];
      if (sign>0.0) {
        for (int i=0; i<3; i++) {
          r0[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i);
          r1[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+3);
          r2[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+6);
          r2_r1[i] = r2[i] - r1[i];
        }
      }
      else {
        for (int i=0; i<3; i++) {
          r0[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+9);
          r1[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+6);
          r2[i] = local_RWGNumbers_trianglesCoord(startIndex_in_localArrays + RWG_index, i+3);
          r2_r1[i] = r2[i] - r1[i];
        }
      }
      for (int i=0; i<3; i++) r_p.push_back(r0[i]);
      const float l_p = sqrt(r2_r1[0]*r2_r1[0] + r2_r1[1]*r2_r1[1] + r2_r1[2]*r2_r1[2]);
      const float RWG_weight = sign * l_p/2.0/sum_weigths;
      TriangleToRWGweightTmp[j][p] = RWG_weight;
    }
    std::vector<float>(r_p).swap(r_p);
    TriangleToRWG_roppTmp.push_back(r_p);
  }
  // now filling the cube values
  TriangleToRWGindex.reserve(N_total_RWG);
  TriangleToRWGweight.reserve(N_total_RWG); 
  TriangleToRWG_ropp.reserve(N_total_RWG*3);
  for (unsigned int j=0; j<T; j++) {
    const int n_rwg = Triangle_numberOfRWGs[j];
    for (int p=0; p<n_rwg; p++) {
      TriangleToRWGindex.push_back(TriangleToRWGindexTmp[j][p]);
      TriangleToRWGweight.push_back(TriangleToRWGweightTmp[j][p]);
      for (int i=0; i<3; i++) TriangleToRWG_ropp.push_back(TriangleToRWG_roppTmp[j][p*3+i]);
    }
  }
}


void Cube::copyCube(const Cube& cubeToCopy) // copy member function
{
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
  for (int i=0; i<3; i++) rCenter[i] = cubeToCopy.rCenter[i];
  for (int i=0; i<3; i++) absoluteCartesianCoord[i] = cubeToCopy.absoluteCartesianCoord[i];
  RWG_numbers.resize(cubeToCopy.RWG_numbers.size());
  RWG_numbers = cubeToCopy.RWG_numbers;
  RWG_numbers_CFIE_OK.resize(cubeToCopy.RWG_numbers_CFIE_OK.size());
  RWG_numbers_CFIE_OK = cubeToCopy.RWG_numbers_CFIE_OK;
  // resize
  Triangle_numberOfRWGs.resize(cubeToCopy.Triangle_numberOfRWGs.size());
  TriangleToRWGindex.resize(cubeToCopy.TriangleToRWGindex.size());
  TriangleToRWGweight.resize(cubeToCopy.TriangleToRWGweight.size());
  TriangleToRWG_ropp.resize(cubeToCopy.TriangleToRWG_ropp.size());
  // copy
  Triangle_numberOfRWGs = cubeToCopy.Triangle_numberOfRWGs;
  TriangleToRWGindex = cubeToCopy.TriangleToRWGindex;
  TriangleToRWGweight = cubeToCopy.TriangleToRWGweight;
  TriangleToRWG_ropp = cubeToCopy.TriangleToRWG_ropp;

  // JPA : shape est plus propre que de le faire dans les 2 dimensions separement
  triangle_GaussCoord.resize(cubeToCopy.triangle_GaussCoord.shape());
  triangle_nHat.resize(cubeToCopy.triangle_nHat.shape());
  triangle_GaussCoord = cubeToCopy.triangle_GaussCoord;
  triangle_nHat = cubeToCopy.triangle_nHat;
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
  Triangle_numberOfRWGs.clear();
  TriangleToRWGindex.clear();
  TriangleToRWGweight.clear();
  TriangleToRWG_ropp.clear();
  triangle_GaussCoord.free();
  triangle_nHat.free();
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

