#ifndef MESH_H
#define MESH_H
#include <blitz/array.h>
#include <vector>
#include <algorithm>

using namespace std;

#include "readWriteBlitzArrayFromFile.h"

// class Mesh {
// 
//   public:
//     int C; // number of cubes
//     int V; // number of vertexes
//     int T; // number of triangles
//     int E; // number of edges
//     int S; // number of surfaces
//     blitz::Array<int, 1> triangles_surfaces;
//     blitz::Array<int, 1> isClosedSurface;
//     blitz::Array<double, 2> cubes_centroids;
//     blitz::Array<double, 2> vertexes_coord;
//     blitz::Array<int, 1> cubes_RWGsNumbers;
//     blitz::Array<int, 2> RWGNumber_signedTriangles;
//     blitz::Array<int, 2> RWGNumber_edgeVertexes;
//     blitz::Array<int, 2> RWGNumber_oppVertexes;
//     blitz::Array<int, 1> RWGNumber_CFIE_OK;
//     // triangles_normals can easily be computed. See vector_functions.cpp
// 
// //  public:
//     // constructors
//     Mesh(void){};
//     Mesh(const string /*path*/);
//     void setMeshFromFile(const string /*path*/); // not really a constructor, but oh well....
//     void copyMesh(const Mesh&);
//     Mesh(const Mesh&); // copy constructor
//     Mesh& operator=(const Mesh&); // copy assignment operator
//     void resizeToZero(void);
//     ~Mesh(void);
// 
//     const int getC(void) const {return C;};
//     void setC(int N) {C = N;};
//     const int getV(void) const {return V;};
//     void setV(int N) {V = N;};
//     const int getT(void) const {return T;};
//     void setT(int N) {T = N;};
//     const int getE(void) const {return E;};
//     void setE(int N) {E = N;};
//     const int getS(void) const {return S;};
//     void setS(int N) {S = N;};
// 
// };

class LocalMesh {

  public:
    int N_local_RWG;
    blitz::Array<int, 1> localRWGNumbers;
    blitz::Array<int, 1> reallyLocalRWGNumbers;
    blitz::Array<int, 1> localRWGNumber_CFIE_OK;
    blitz::Array<float, 2> localRWGNumber_trianglesCoord;

    // constructors
    LocalMesh(void);
    LocalMesh(const string /*path*/);
    void setLocalMeshFromFile(const string /*path*/); // not really a constructor, but oh well....
    void copyLocalMesh(const LocalMesh&);
    LocalMesh(const LocalMesh&); // copy constructor
    LocalMesh& operator=(const LocalMesh&); // copy assignment operator
    ~LocalMesh(void);
    void resizeToZero(void);
    void writeLocalMeshToFile(const string /*path*/);
};

// we need the following Dictionary class in various places
template <typename K, typename V>
class Dictionary {
    K key;
    V val;
  public:
    // constructors
    Dictionary(const K k, const V v) {key = k; val = v;};
    const V getVal(void) const {return val;};
    const K getKey(void) const {return key;};
    Dictionary(const Dictionary & dictionaryToCopy) {val = dictionaryToCopy.getVal(); key = dictionaryToCopy.getKey();}; // copy constructor
    ~Dictionary(){};
    // overloaded operators for sorting
    bool operator== (const Dictionary & right) const {if ( getKey() == right.getKey() ) return 1; else return 0;};
    bool operator< (const Dictionary & right) const {if ( getKey() < right.getKey() ) return 1; else return 0;};
};


void compute_RWGNumber_trianglesNumbers(blitz::Array<int, 1>& indexes_of_triangles,
                                        const blitz::Array<int, 1>& list_of_RWG_numbers,
                                        const blitz::Array<int, 2>& RWGNumber_signedTriangles);

void compute_indexesEqualEdges(std::vector<std::vector<int> >& indexesEqualEdges,
                               const blitz::Array<int, 1>& indexesEqualPreceding,
                               const blitz::Array<int, 1>& ind_sorted_e_v);

void compute_edgeNumber_vertexes(blitz::Array<int, 2>& edgeNumber_vertexes,
                                 const std::vector<std::vector<int> >& indexesEqualEdges,
                                 const blitz::Array<int, 2>& col_sorted_e_v);

void compute_edgeNumber_triangles(blitz::Array<int, 2>& edgeNumber_triangles,
                                  const std::vector<std::vector<int> >& indexesEqualEdges);

void compute_triangle_adjacentTriangles(std::vector<std::vector<int> >& triangle_adjacentTriangles,
                                        const blitz::Array<int, 2>& edgeNumber_triangles,
                                        const int T);

void reorder_triangle_vertexes(blitz::Array<int, 2>& triangle_vertexes,
                               blitz::Array<int, 1>& triangles_surfaces,
                               const blitz::Array<double, 2>& vertexes_coord,
                               const std::vector<std::vector<int> >& triangle_adjacentTriangles);

void is_surface_closed(blitz::Array<int, 1>& is_closed_surface,
                       blitz::Array<std::vector<int>, 2>& connected_surfaces,
                       blitz::Array<std::vector<int>, 2>& potential_closed_surfaces,
                       const blitz::Array<int, 1>& triangles_surfaces,
                       const blitz::Array<int, 2>& edgeNumber_triangles);

void RWGNumber_signedTriangles_computation(blitz::Array<int, 2>& RWGNumber_signedTriangles,
                                           blitz::Array<int, 2>& RWGNumber_edgeVertexes,
                                           const blitz::Array<int, 2>& edgeNumber_triangles,
                                           const blitz::Array<int, 2>& edgeNumber_vertexes,
                                           const blitz::Array<int, 1>& triangles_surfaces,
                                           const blitz::Array<int, 1>& is_closed_surface,
                                           const blitz::Array<int, 2>& triangle_vertexes,
                                           const blitz::Array<double, 2>& vertexes_coord);

void RWGNumber_oppVertexes_computation(blitz::Array<int, 2>& RWGNumber_oppVertexes,
                                       const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                                       const blitz::Array<int, 2>& RWGNumber_edgeVertexes,
                                       const blitz::Array<int, 2>& triangle_vertexes);


#endif
