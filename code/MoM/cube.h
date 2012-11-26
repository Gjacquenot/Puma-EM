#ifndef CUBE_H
#define CUBE_H

#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <vector>
#include <algorithm>  // Include STL algorithms for sorting lists and vectors

using namespace std;

#include "FMM.h"

class Cube {
    //! tells if the current cube is leaf (finest) or not
    bool leaf;
    //! the absolute number of the current cube
    int number;
    //! the current index of the current cube in the list of all cubes at a given level
    int index;
    //! index before sorting
    int oldIndex;
    //! the absolute number of the father
    int fatherNumber;
    //! index of the father 
    int fatherIndex;
  public:
    //! the process number of the current cube
    int procNumber;
    //! the process number of the father of the current cube
    int fatherProcNumber;
    //! the indexes of the children of the current cube
    std::vector<int> sonsIndexes;
    //! the processes number of the children of the current cube
    std::vector<int> sonsProcNumbers;
    //! the indexes of the neighbors of the current cube
    std::vector<int> neighborsIndexes;
    //! \brief the indexes of the cubes that will participate to the alpha translation to the current cube
    //! and that are located on the same process as that of the current cube
    std::vector<int> localAlphaTransParticipantsIndexes;
    //! \brief the indexes of the cubes that will participate to the alpha translation to the current cube
    //! and that are located on another process as that of the current cube
    std::vector<int> nonLocalAlphaTransParticipantsIndexes;
    //! the absolute coordinates of the center of the current cube 
    float rCenter[3];
    //! \brief the absolute cartesian coordinates of the center of the current cube
    //! calculated with respect to the location of the Father of all Cubes and the side length of the current cube  
    float absoluteCartesianCoord[3];
    //! the numbers of the RWGs pertaining to the current cube
    std::vector<int> RWG_numbers;
    std::vector<int> RWG_numbers_CFIE_OK;
    blitz::Array<blitz::TinyVector<float, 3>, 2> GaussLocatedWeightedRWG;
    blitz::Array<blitz::TinyVector<float, 3>, 2> GaussLocatedWeighted_nHat_X_RWG;
    blitz::Array<blitz::TinyVector<float, 3>, 2> GaussLocatedExpArg;

    // constructors
    Cube(void){};
    //! the base constructor
    /*!
      \param is_leaf a bool argument, true if cube is leaf
      \param level the number of the level
      \param sideLength the length of the cube side
      \param big_cube_lower_coord the lower coordinate of the father of all cubes
      \param r_c the coordinates of cube center
      \return the constructed cube
    */
    Cube(const bool is_leaf, 
         const int level,
         const double sideLength,
         const double big_cube_lower_coord[3],
         const blitz::Array<double, 1>& r_c);
    //! allows a cube to be overwritten by the provided cube
    /*!
      \param cubeToCopy the cube to copy
    */
    void copyCube (const Cube& cubeToCopy);
    Cube(const Cube&); // copy constructor
    Cube& operator=(const Cube&); // copy assignment operator
    //! the constructor from a children cube
    /*!
      \param sonCube the children cube
      \param level the number of the level
      \param big_cube_lower_coord the lower coordinate of the father of all cubes
      \param sideLength the length of the cube side
      \return the constructed cube
    */
    Cube(const Cube& sonCube,
         const int level,
         const double big_cube_lower_coord[3],
         const double sideLength);
    //! the destructor
    ~Cube();

    // specific functions
    //! returns true if we are at leaf level, false otherwise
    const bool getLeaf(void) const {return leaf;}
    void addSon(const Cube&);
    //! returns the size of the sonsIndexes vector
    const int getSonsIndexesSize(void) const {return sonsIndexes.size();};
    //! returns the sonsIndexes vector
    const vector<int> getSonsIndexes(void) const {return sonsIndexes;};
    //! returns the porocesses numbers of the sons, arranged as a vector
    const vector<int> getSonsProcNumbers(void) const {return sonsProcNumbers;};
    const vector<int> getNeighborsIndexes(void) const {return neighborsIndexes;};
    const vector<int> getLocalAlphaTransParticipantsIndexes(void) const {return localAlphaTransParticipantsIndexes;};
    const vector<int> getNonLocalAlphaTransParticipantsIndexes(void) const {return nonLocalAlphaTransParticipantsIndexes;};
//    const blitz::TinyVector<float, 3> getRCenter(void) const {return rCenter;};
//    const blitz::TinyVector<float, 3> getAbsoluteCartesianCoord(void) const {return absoluteCartesianCoord;};
    const vector<int> getRWG_numbers(void) const {return RWG_numbers;};
    const vector<int> getRWG_numbers_CFIE_OK(void) const {return RWG_numbers_CFIE_OK;};
    void setIndex(const int i) {index = i;};
    const int getIndex(void) const {return index;};
    void setOldIndex(const int i) {oldIndex = i;};
    const int getOldIndex(void) const {return oldIndex;};
    const int getNumber(void) const {return number;};
    const int getProcNumber(void) const {return procNumber;};
    const int getFatherNumber(void) const {return fatherNumber;};
    const int getFatherProcNumber(void) const {return fatherProcNumber;};
    void setFatherNumber(const int n) {fatherNumber = n;};
    const int getFatherIndex(void) const {return fatherIndex;};
    void setFatherIndex(const int i) {fatherIndex = i;};
    const int getSonIndex(const int i) {return sonsIndexes[i];};

    //! \brief computes the points locations and values for the arguments for the complex exponentials 
    //! that will be used in computing the radiation function of the leaf cube 
    /*!
      \param target_mesh the mesh of the target, C++ object
      \param N_Gauss int, the number of points per triangles
      \return void
    */
    void computeGaussLocatedArguments(const blitz::Array<int, 1>& local_RWG_numbers,
                                      const blitz::Array<int, 1>& local_RWG_Numbers_CFIE_OK,
                                      const blitz::Array<float, 2>& local_RWGNumbers_trianglesCoord,
                                      const int startIndex_in_localArrays,
                                      const int NRWG,
                                      const int N_Gauss);
    const blitz::Array<blitz::TinyVector<float, 3>, 2> getGaussLocatedWeightedRWG(void) const {return GaussLocatedWeightedRWG;};
    const blitz::Array<blitz::TinyVector<float, 3>, 2> getGaussLocatedWeighted_nHat_X_RWG(void) const {return GaussLocatedWeighted_nHat_X_RWG;};
    const blitz::Array<blitz::TinyVector<float, 3>, 2> getGaussLocatedExpArg(void) const {return GaussLocatedExpArg;};

    // overloaded operators
    bool operator== (const Cube &) const;
    bool operator< (const Cube &) const;
};


#endif
