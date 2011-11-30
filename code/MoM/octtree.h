#ifndef OCTTREE_H
#define OCTTREE_H

#include <iostream>
#include <string>
#include <complex>
#include <blitz/array.h>
#include <vector>
#include <algorithm>  // Include STL algorithms for sorting lists and vectors

using namespace std;

#include "mesh.h"
#include "level.h"

/*! \class Octtree
    \brief the Octtree class is basically a vector container holding the levels, and additional information and arrays.

    The Octtree class is composed of:
      - <strong> vector<Cube> levels </strong>: a vector container holding all the levels of the tree
      - ...
*/

class Octtree {
    int L;
    int numberOfUpdates;
    int procNumber;
    int totalNumProcs;
  public:
    complex<double> k;
    complex<float> eps_r;
    complex<float> mu_r;
    float w;
    int VERBOSE;
    int N_GaussOnTriangle;
    int DIRECTIONS_PARALLELIZATION, N_levels, ALLOW_CEILING_LEVEL;
    vector<Level> levels;
    blitz::TinyVector<double, 3> big_cube_lower_coord;
    blitz::TinyVector<double, 3> big_cube_center_coord;
    blitz::Array<std::complex<float>, 1> CFIE;
    string octtreeDataPath;

    // constructors
    Octtree(void) {};
    Octtree(const string /*octtree_data_path*/,
            const blitz::Array<double, 2>& /*cubes_centroids*/,
            const int /*proc_id*/,
            const int /*num_procs*/);
    void copyOcttree (const Octtree&);
    Octtree(const Octtree &); // copy constructor
    Octtree& operator=(const Octtree&); // copy assignment operator
    void constructArrays(void);
    void computeGaussLocatedArguments(const blitz::Array<int, 1>& /*local_cubes_NRWG*/, 
                                      const blitz::Array<int, 1>& /*local_RWG_numbers*/, 
                                      const blitz::Array<int, 1>& /*local_RWG_Numbers_CFIE_OK*/, 
                                      const blitz::Array<float, 2>& /*local_RWGNumbers_trianglesCoord*/);
    void RWGs_renumbering(void);
    void computeIndexesOfCubesInOriginalMesh(blitz::Array<int, 1>& /*oldIndexesOfCubes*/);
    ~Octtree();

    // functions
    const int getNumberOfUpdates (void) const {return numberOfUpdates;};
    const int setNumberOfUpdates (const int n) {numberOfUpdates = n;};
    const int getL(void) const {return L;};
    void setProcNumber (const int n) {procNumber = n;};
    const int getProcNumber (void) const {return procNumber;};
    void setTotalNumProcs (const int n) {totalNumProcs = n;};
    const int getTotalNumProcs (void) const {return totalNumProcs;};
    const int getLevelsSize (void) const {return levels.size();};
    const complex<double> getK(void) const {return k;};
    const complex<float> getEps_r(void) const {return eps_r;};
    const complex<float> getMu_r(void) const {return mu_r;};
    const float getW(void) const {return w;};
    const Level getLevel(const int l) const {return levels[l];};
    const Cube getCubeLevel(const int i, const int l) const {return levels[l].getCube(i);}; 
    void assignCubesToProcessors(const int /*num_procs*/, const int /*CUBES_DISTRIBUTION*/);
    void writeAssignedLeafCubesToDisk(const string /*path*/, const string /*filename*/);
    void updateSup(const blitz::Array<std::complex<float>, 1>&); // coefficients of RWG functions
    void exchangeSupsIndividually(blitz::Array< blitz::Array<std::complex<float>, 2>, 1>& /*SupThisLevel*/, const int /*l*/, const vector<int> & /*localCubesIndexes*/);
    void exchangeSupsInBlocks(blitz::Array< blitz::Array<std::complex<float>, 2>, 1>& /*SupThisLevel*/, const int /*l*/, const vector<int> & /*localCubesIndexes*/);
    void ZIFarComputation(blitz::Array<std::complex<float>, 1>& /*ZI*/, // result of matrix-vector multiplication
                          const blitz::Array<std::complex<float>, 1>& /*I_PQ*/); // coefficients of RWGs
    const blitz::Array<std::complex<float>, 1> getCFIE(void) const {return CFIE;};
    const blitz::Array<complex<float>, 1> getAlphaLevel(const int m, const int n, const int p, const int l) { return levels[l].getAlpha( m, n, p);};
    const vector<int> getNeighborsSonsIndexes(const int, const int) const;
    void findAlphaTransParticipantsIndexes(const int l);
    void SupAlphaMultiplication(blitz::Array<std::complex<float>, 2>& /*SupAlpha*/,
                                const blitz::Array<std::complex<float>, 2>& /*Sup*/,
                                const blitz::Array<std::complex<float>, 1>& /*alphaTranslation*/,
                                const blitz::Array<int, 1>& /*alphaTranslationIndexesNonZeros*/,
                                const blitz::Array<int, 1>& /*alphaTranslationIndexes*/,
                                const blitz::TinyVector<int, 3>& /*alphaCartesianCoord*/);
    void SupAlphaMultiplicationDirections(blitz::Array<std::complex<float>, 2>& /*SupAlpha*/,
                                                   const blitz::Array<std::complex<float>, 2>& /*Sup*/,
                                                   const blitz::Array<std::complex<float>, 1>& /*alphaTranslation*/,
                                                   const blitz::Array<int, 1>& /*alphaTranslationIndexesNonZeros*/,
                                                   const blitz::TinyVector<int, 3>& /*alphaCartesianCoord*/);
    void alphaTranslationsToCube(blitz::Array<std::complex<float>, 2>& /*S_tmp*/,
                                 const blitz::Array< blitz::Array<std::complex<float>, 2>, 1>& /*LevelSup*/,
                                 const int /*l*/,
                                 const int /*cubeIndex*/,
                                 const std::vector<int>& /*indexesAlphaParticipants*/,
                                 const int /*DIRECTIONS_PARALLELIZATION*/);
    void shiftExp(blitz::Array<std::complex<float>, 2> /*S*/,
                  const blitz::Array<std::complex<float>, 1>& /*shiftingArray*/);
    void S2DWeighting (blitz::Array<std::complex<float>, 2> /*S*/,
                       const blitz::Array<float, 1>& /*Wtheta*/,
                       const blitz::Array<float, 1>& /*Wphi*/);
    void computeFarField (blitz::Array<std::complex<float>, 2>& e_theta_far,
                          blitz::Array<std::complex<float>, 2>& e_phi_far,
                          const blitz::Array<float, 1>& octtreeXthetas_coarsest,
                          const blitz::Array<float, 1>& octtreeXphis_coarsest,
                          const blitz::Array<std::complex<float>, 1>& I_PQ,
                          const string octtree_data_path);
    void resizeSdownLevelsToZero(void) {for (int i=0 ; i<levels.size() ; ++i) levels[i].Sdown.resize(0);};
};

#endif
