#ifndef LEVEL_H
#define LEVEL_H

#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <vector>
#include <algorithm>  // Include STL algorithms for sorting lists and vectors

using namespace std;

#include "dictionary.h"
#include "cube.h"
#include "interpolation.h"

/*! \class Level
    \brief the Level class is basically a vector container holding the cubes pertaining to this level, and additional information and arrays.

    The Level class is composed of:
      - <strong> vector<Cube> cubes </strong>: a vector container holding all the cubes of the level
      - <strong> vector<Dictionary> numbersToIndexes </strong>:  object that returns the index for a given cube number
      - <strong> vector<int> localCubesIndexes </strong>: the list of the cubes local to the current process
      - <strong> vector< vector<int> > listOfFcToBeReceived </strong>: a list of radiation functions to be received from each process
      - <strong> vector< vector<int> > listOfFcToBeSent </strong>: a list of radiation functions to be sent to each process
      - <strong> vector<int> cubesIndexesAfterReduction </strong>: the indexes of the cubes after removal of the nonlocal cubes 
*/
class Level {
    //! tells if the current level is leaf (finest) or not
    bool leaf;
    //! tells if the current level is ceiling (the last one)
    bool ceiling;
    //! the level number
    int level;
    //! the maximum number of subdivisions of the father of all cubes 
    int maxNumberCubes1D;
    //! the numbers of cubes in direction X
    int NCubesX;
    //! the numbers of cubes in direction Y
    int NCubesY;
    //! the numbers of cubes in direction Z
    int NCubesZ; 
    //! the offsets for alpha indexes
    int offsetAlphaIndexX, offsetAlphaIndexY, offsetAlphaIndexZ;
    //! the number of terms in the alpha expansion
    int N;
    //! the length of the side of the cubes that belong to the current level
    double cubeSideLength;
    //! the wavenumber
    std::complex<double> k;
    //! the number of times it has been copied. Important data/statistics for code optimization reasons
    int numberTimesCopied;
  public:
    //! a list of the cubes that belong to the current level 
    std::vector<Cube> cubes;
    //! a cube cartesian number to cube index array. Necessary to quickly find its neighbours
    std::vector< Dictionary<int, int> > numbersToIndexes;
    //! a list of the indexes of the cubes that are on the same process 
    std::vector<int> localCubesIndexes;
    //! the list of indexes of the cubes whose radiation functions must be received
    /*!
      "listOfFcToBeSent" and "listOfFcToBeReceived". Each list is as long as the number of processes, 
      and each element is accessed by the process number with which the communication will be 
      established. Each element consists of a list of indexes corresponding to the radiation 
      functions to be sent or received by the current process at the current level.
    */
    std::vector< std::vector<int> > listOfFcToBeReceived;
    //! the list of indexes of the cubes whose radiation functions must be sent
    std::vector< std::vector<int> > listOfFcToBeSent;
    //! the indexes of the cubes that remain after the level has been resized to hold only the local cubes
    std::vector<int> cubesIndexesAfterReduction;
    //! the theta sampling of the current level 
    blitz::Array<float, 1> thetas;
    //! the phi sampling of the current level 
    blitz::Array<float, 1> phis;
    //! the theta weights of the current level 
    blitz::Array<float, 1> weightsThetas;
    //! the phi weights of the current level 
    blitz::Array<float, 1> weightsPhis;
    blitz::Array<std::complex<float>, 2> shiftingArrays;
    blitz::Array< blitz::Array<std::complex<float>, 1>, 3> alphaTranslations;
    blitz::Array< blitz::Array<int, 1>, 3> alphaTranslationsIndexesNonZeros;
    blitz::Array<int, 4> alphaTranslationsIndexes; // necessary due to the use of symmetry
    LagrangeFastInterpolator2D lfi2D; // the interpolator for the next level
    blitz::Array< blitz::Array<std::complex<float>, 2>, 1> Sdown;
    //! tells if we have parallelization by directions (currently only for the ceiling level)
    int DIRECTIONS_PARALLELIZATION;
    blitz::Array<int, 1> MPI_Scatterv_scounts, MPI_Scatterv_displs;

//  public:
    // constructors
    Level();
    //! a simplified constructor for the leaf level
    /*!
      \param l the level number
      \param leaf_side_length the side length of the leaf cubes
      \param big_cube_lower_coord the lower coordinate of the father of all cubes
      \param cubes_centroids the coordinates of the centroids for each cube
      \return the constructed level
    */
    Level(const int l,
          const double leaf_side_length,
          const double big_cube_lower_coord[3],
          const blitz::Array<double, 2>& cubes_centroids);
    //! the complete constructor for the leaf level
    /*!
      \param l the level number
      \param N_expansion the number of terms in the expansion of the Green's function
      \param leaf_side_length the side length of the leaf cubes
      \param big_cube_lower_coord the lower coordinate of the father of all cubes
      \param cubes_centroids the coordinates of the centroids for each cube
      \param wavenumber the wavenumber
      \param A_theta the lower limit for \f$ \theta \f$ interval : \f$ 0 \f$ 
      \param B_theta the upper limit for \f$ \theta \f$ interval : \f$ \pi \f$ 
      \param Xthetas the \f$ \theta \f$  sampling abscissas
      \param Wthetas the \f$ \theta \f$  sampling weights
      \param INCLUDED_THETA_BOUNDARIES 1 if \f$ \theta \f$ interval boundaries are included in the sampling, 0 otherwise
      \param N_theta number of \f$ \theta \f$ samples
      \param PERIODIC_Theta 1 if we use Gaussian interpolator, 0 if we want to use the Lagrangian interpolator
      \param CYCLIC_Theta 1 if we use cyclic interpolation for the radiation functions, 0 otherwise
      \param NOrderInterpolatorTheta the order of the interpolator
      \param A_phi the lower limit for \f$ \phi \f$ interval : \f$ 0 \f$ 
      \param B_phi the upper limit for \f$ \phi \f$ interval : \f$ 2 \pi \f$ 
      \param Xphis the \f$ \phi \f$ sampling abscissas
      \param Wphis the \f$ \phi \f$ sampling weights
      \param INCLUDED_PHI_BOUNDARIES 1 if \f$ \phi \f$ boundaries are included in the sampling, 0 otherwise
      \param N_phi number of \f$ \phi \f$ samples
      \param PERIODIC_Phi 1 if we use Gaussian interpolator, 0 if we want to use the Lagrangian interpolator
      \param CYCLIC_Phi 1 if we use cyclic interpolation for the radiation functions, 0 otherwise
      \param NOrderInterpolatorPhi the order of the interpolator
      \param XthetasNextLevel the \f$ \theta \f$ sampling abscissas of the parent level
      \param XphisNextLevel the \f$ \phi \f$ sampling abscissas of the parent level
      \param VERBOSE if 1, detailed output

      \return the constructed level
    */
    Level(const int l,
          const int N_expansion,
          const double leaf_side_length,
          const double big_cube_lower_coord[3],
          const blitz::Array<double, 2>& cubes_centroids,
          const std::complex<double>& wavenumber,
          const float A_theta,
          const float B_theta,
          const blitz::Array<float, 1>& Xthetas,
          const blitz::Array<float, 1>& Wthetas,
          const int INCLUDED_THETA_BOUNDARIES,
          const int N_theta,
          const int PERIODIC_Theta,
          const int CYCLIC_Theta,
          const int NOrderInterpolatorTheta,
          const float A_phi,
          const float B_phi,
          const blitz::Array<float, 1>& Xphis,
          const blitz::Array<float, 1>& Wphis,
          const int INCLUDED_PHI_BOUNDARIES,
          const int N_phi,
          const int PERIODIC_Phi,
          const int CYCLIC_Phi,
          const int NOrderInterpolatorPhi,
          const blitz::Array<float, 1>& XthetasNextLevel,
          const blitz::Array<float, 1>& XphisNextLevel,
          const int VERBOSE);
    //! allows a level to be overwritten by the provided level
    /*!
      \param levelToCopy the level to copy
    */
    void copyLevel (const Level& levelToCopy);
    Level(const Level &); // copy constructor
    Level& operator=(const Level&); // copy assignment operator
    //! the constructor from a child level
    /*!
      \param sonLevel the child level
      \param N_expansion the number of terms in the expansion of the Green's function
      \param big_cube_lower_coord the lower coordinate of the father of all cubes
      \param A_theta the lower limit for \f$ \theta \f$ interval : \f$ 0 \f$ 
      \param B_theta the upper limit for \f$ \theta \f$ interval : \f$ \pi \f$ 
      \param Xthetas the \f$ \theta \f$  sampling abscissas
      \param Wthetas the \f$ \theta \f$  sampling weights
      \param INCLUDED_THETA_BOUNDARIES 1 if \f$ \theta \f$ interval boundaries are included in the sampling, 0 otherwise
      \param N_theta number of \f$ \theta \f$ samples
      \param PERIODIC_Theta 1 if we use Gaussian interpolator, 0 if we want to use the Lagrangian interpolator
      \param CYCLIC_Theta 1 if we use cyclic interpolation for the radiation functions, 0 otherwise
      \param NOrderInterpolatorTheta the order of the interpolator
      \param A_phi the lower limit for \f$ \phi \f$ interval : \f$ 0 \f$ 
      \param B_phi the upper limit for \f$ \phi \f$ interval : \f$ 2 \pi \f$ 
      \param Xphis the \f$ \phi \f$ sampling abscissas
      \param Wphis the \f$ \phi \f$ sampling weights
      \param INCLUDED_PHI_BOUNDARIES 1 if \f$ \phi \f$ boundaries are included in the sampling, 0 otherwise
      \param N_phi number of \f$ \phi \f$ samples
      \param PERIODIC_Phi 1 if we use Gaussian interpolator, 0 if we want to use the Lagrangian interpolator
      \param CYCLIC_Phi 1 if we use cyclic interpolation for the radiation functions, 0 otherwise
      \param NOrderInterpolatorPhi the order of the interpolator
      \param XthetasNextLevel the \f$ \theta \f$ sampling abscissas of the parent level
      \param XphisNextLevel the \f$ \phi \f$ sampling abscissas of the parent level
      \param VERBOSE if 1, detailed output

      \return the constructed level
    */
    Level(const Level & sonLevel,
          const int N_expansion,
          const double big_cube_lower_coord[3],
          const float A_theta,
          const float B_theta,
          const blitz::Array<float, 1>& Xthetas,
          const blitz::Array<float, 1>& Wthetas,
          const int INCLUDED_THETA_BOUNDARIES,
          const int N_theta,
          const int PERIODIC_Theta,
          const int CYCLIC_Theta,
          const int NOrderInterpolatorTheta,
          const float A_phi,
          const float B_phi,
          const blitz::Array<float, 1>& Xphis,
          const blitz::Array<float, 1>& Wphis,
          const int INCLUDED_PHI_BOUNDARIES,
          const int N_phi,
          const int PERIODIC_Phi,
          const int CYCLIC_Phi,
          const int NOrderInterpolatorPhi,
          const blitz::Array<float, 1>& XthetasNextLevel,
          const blitz::Array<float, 1>& XphisNextLevel,
          const int VERBOSE); // fromSon constructor
    //! the destructor
    ~Level();

    // specific functions
    bool getLeaf(void) const {return leaf;}
    void setCeiling(bool l) {ceiling = l;}
    bool getCeiling(void) const {return ceiling;}
    void setLevel(int l) {level = l;}
    int getLevel(void) const {return level;}
    //! maxNumberCubes1D = static_cast<int> pow(2.0, level)
    int getMaxNumberCubes1D(void) const {return maxNumberCubes1D;}
    int getNumberTimesCopied(void) const {return numberTimesCopied;}
    void incrementNumberTimesCopied(void) {numberTimesCopied += 1;}
    void NCubesXYZComputation(const int VERBOSE);
    int getNCubesX(void) const {return NCubesX;}
    int getNCubesY(void) const {return NCubesY;}
    int getNCubesZ(void) const {return NCubesZ;}
    int getOffsetAlphaIndexX(void) const {return offsetAlphaIndexX;}
    int getOffsetAlphaIndexY(void) const {return offsetAlphaIndexY;}
    int getOffsetAlphaIndexZ(void) const {return offsetAlphaIndexZ;}
    int getN(void) const {return N;}
    double getCubeSideLength(void) const {return cubeSideLength;}
    std::complex<double> getK(void) const {return k;}
    void addNode(Cube cube) {cubes.push_back(cube);}
    Cube getCube(const int i) const {return cubes[i];}
    std::vector<Cube> getCubes(void) const {return cubes;}
    double getCubesSizeMB(void) const {return cubes.size() * (thetas.size()*phis.size()) *  2 * 2.0*4.0/(1024.0*1024.0);}
    std::vector<int> getCubeNeighbors(const int i) const {return cubes[i].getNeighborsIndexes();}
    std::vector< Dictionary<int, int> > getNumbersToIndexes(void) const {return numbersToIndexes;}
    std::vector<int> getLocalCubesIndexes(void) const {return localCubesIndexes;}
    std::vector< std::vector<int> > getListOfFcToBeReceived(void) const {return listOfFcToBeReceived;}
    std::vector< std::vector<int> > getListOfFcToBeSent(void) const {return listOfFcToBeSent;}
    std::vector<int> getCubesIndexesAfterReduction(void) const {return cubesIndexesAfterReduction;}
    void computeOldIndexesOfCubes(blitz::Array<int, 1>& /*oldIndexesOfCubes*/); 
    void computeLevelReduction(void);
    int getNumbersToIndexesSize(void) const {return numbersToIndexes.size();}
    int getNumberToIndex(const int i) const {return numbersToIndexes[i].getKey();}
    int getIndexToIndex(const int i) const {return numbersToIndexes[i].getVal();}
    int getIndexOfNumber(const int) const;
    int getLevelSize(void) const {return cubes.size();}
    int getSizeOfAlphaTransParticipantsIndexes(void) const {int result = 0; for (unsigned int i=0; i<cubes.size(); ++i) result += (cubes[i].localAlphaTransParticipantsIndexes.size() + cubes[i].nonLocalAlphaTransParticipantsIndexes.size()); return result;}
    float getSizeMBOfAlphaTransParticipantsIndexes(void) const {return getSizeOfAlphaTransParticipantsIndexes()*4.0/(1024.0*1024.0);}
    blitz::Array<float, 1> getThetas(void) const {return thetas;}
    blitz::Array<float, 1> getPhis(void) const {return phis;}
    int getNThetas(void) const {return thetas.size();}
    int getNPhis(void) const {return phis.size();}
    void computeGaussLocatedArguments(const blitz::Array<int, 1>& /*local_cubes_NRWG*/, 
                                      const blitz::Array<int, 1>& /*local_RWG_numbers*/, 
                                      const blitz::Array<int, 1>& /*local_RWG_Numbers_CFIE_OK*/, 
                                      const blitz::Array<int, 2>& /*local_RWGNumbers_signedTriangles*/,
                                      const blitz::Array<float, 2>& /*local_RWGNumbers_trianglesCoord*/,
                                      const int /*N_Gauss*/);
    void RWGs_renumbering(void);
    void shiftingArraysComputation(void);
    double getShiftingArraysSizeMB(void) const {return shiftingArrays.size() *  2.0*4.0/(1024.0*1024.0);}
    blitz::Array<int, 1> getAlphaTranslationsExtents(void) const;
    blitz::Array< blitz::Array<std::complex<float>, 1>, 3> getAlphaTranslations(void) const {return alphaTranslations;}
    double getAlphaTranslationsSizeMB(void) const;
    double getAlphaTranslationsIndexesSizeMB(void) const {return alphaTranslationsIndexes.size() * 4.0/(1024.0*1024.0);}
    blitz::Array<std::complex<float>, 2> getShiftingArrays(void) const {return shiftingArrays;}
    const blitz::Array<std::complex<float>, 1> getShiftingArray(const double Dx, const double Dy, const double Dz) const {return shiftingArrays((Dx>0.0) * 4 + (Dy>0.0) * 2 + (Dz>0.0), blitz::Range::all());}
    void alphaTranslationsComputation(const int VERBOSE,
                                      const float alphaTranslation_smoothing_factor,
                                      const float alphaTranslation_thresholdRelValueMax,
                                      const float alphaTranslation_RelativeCountAboveThreshold);
    void alphaTranslationIndexConstructionZ(blitz::Array<int, 1>& newAlphaIndex,
                                            const blitz::Array<int, 1>& oldAlphaIndex,
                                            const int alphaCartesianCoordZ,
                                            const int N_theta,
                                            const int N_phi);
    void alphaTranslationIndexConstructionY(blitz::Array<int, 1>& /*newAlphaIndex*/,
                                            const blitz::Array<int, 1>& /*oldAlphaIndex*/,
                                            const int /*alphaCartesianCoordY*/,
                                            const int /*N_theta*/,
                                            const int /*N_phi*/);
    void alphaTranslationIndexConstructionX(blitz::Array<int, 1>& /*newAlphaIndex*/,
                                            const blitz::Array<int, 1>& /*oldAlphaIndex*/,
                                            const int /*alphaCartesianCoordX*/,
                                            const int /*N_theta*/,
                                            const int /*N_phi*/);
    blitz::Array<float, 1> getWeightsThetas(void) const {return weightsThetas;} 
    blitz::Array<float, 1> getWeightsPhis(void) const {return weightsPhis;} 
    LagrangeFastInterpolator2D getLfi2D(void) const {return lfi2D;}

    void sortCubesByParents(void);
    void printCubesSonsIndexes(void);
    void printCubesFathersNumbers(void);
    void printCubesRCenters(void);
    void printCubesRWG_numbers(void);
    void printThetaPhi(void);
    void updateFatherIndexes(const Level& fatherLevel);
    void setNumbersToIndexes(void);
    void sortNumbersToIndexes(void);
    void printNumbersToIndexes(void);
    void searchCubesNeighborsIndexes(void);
    void printCubesNeighborsIndexes(void);
    void computeLocalCubesIndexes(const int /*procNumber*/);

    void computeSup(blitz::Array<std::complex<float>, 2> & /*Sup*/,
                    const std::complex<double>& /*k*/,
                    const blitz::Array<std::complex<float>, 1>& /*I_PQ*/,
                    const Cube & /*cube*/,
                    const blitz::Array<float, 1>& /*thetas*/,
                    const blitz::Array<float, 1>& /*phis*/);
    void sphericalIntegration(blitz::Array<std::complex<float>, 1>& /*ZI*/,
                              const blitz::Array<std::complex<float>, 2>& /*Sdown*/,
                              const Cube & /*cube*/,
                              const blitz::Array<float, 1>& /*thetas*/,
                              const blitz::Array<float, 1>& /*phis*/,
                              const float /*w*/,
                              const std::complex<float>& /*mu_r*/,
                              const std::complex<double>& /*k*/,
                              const blitz::Array<std::complex<float>, 1>& /*CFIE*/);
};

#endif
