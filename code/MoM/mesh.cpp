#include "mesh.h"
#include "dictionary.h"

/*Mesh::Mesh(const string path) {setMeshFromFile(path);};

void Mesh::copyMesh(const Mesh& meshToCopy) /// copy member function
{
  C = meshToCopy.C;
  V = meshToCopy.V;
  T = meshToCopy.T;
  E = meshToCopy.E;
  S = meshToCopy.S;
  vertexes_coord.resize(V, 3);
  cubes_centroids.resize(C, 3);
  cubes_RWGsNumbers.resize(meshToCopy.cubes_RWGsNumbers.size());
  RWGNumber_signedTriangles.resize(meshToCopy.RWGNumber_signedTriangles.extent(0), meshToCopy.RWGNumber_signedTriangles.extent(1));
  RWGNumber_edgeVertexes.resize(meshToCopy.RWGNumber_edgeVertexes.extent(0), meshToCopy.RWGNumber_edgeVertexes.extent(1));
  RWGNumber_oppVertexes.resize(meshToCopy.RWGNumber_oppVertexes.extent(0), meshToCopy.RWGNumber_oppVertexes.extent(1));
  RWGNumber_CFIE_OK.resize(meshToCopy.RWGNumber_CFIE_OK.size());
  triangles_surfaces.resize(T);
  isClosedSurface.resize(S);
  vertexes_coord = meshToCopy.vertexes_coord;
  cubes_centroids = meshToCopy.cubes_centroids;
  cubes_RWGsNumbers = meshToCopy.cubes_RWGsNumbers;
  triangles_surfaces = meshToCopy.triangles_surfaces;
  isClosedSurface = meshToCopy.isClosedSurface;
  RWGNumber_signedTriangles = meshToCopy.RWGNumber_signedTriangles;
  RWGNumber_edgeVertexes = meshToCopy.RWGNumber_edgeVertexes;
  RWGNumber_oppVertexes = meshToCopy.RWGNumber_oppVertexes;
  RWGNumber_CFIE_OK = meshToCopy.RWGNumber_CFIE_OK;
}

Mesh::Mesh(const Mesh& meshToCopy) /// copy constructor
{
  copyMesh(meshToCopy);
}

Mesh& Mesh::operator=(const Mesh& meshToCopy) { /// copy assignment
  copyMesh(meshToCopy);
  return *this;
}

Mesh::~Mesh() {
  //cout << "Deleting mesh object..." << endl;
  resizeToZero();
  vertexes_coord.free();
  cubes_centroids.free();
  cubes_RWGsNumbers.free();
  triangles_surfaces.free();
  isClosedSurface.free();
  RWGNumber_signedTriangles.free();
  RWGNumber_edgeVertexes.free();
  RWGNumber_oppVertexes.free();
  RWGNumber_CFIE_OK.free();
}

void Mesh::setMeshFromFile(const string path)
{
  {
    string filename = path + "C.txt";
    readIntFromASCIIFile(filename, C);
  }
  {
    string filename = path + "V.txt";
    readIntFromASCIIFile(filename, V);
  }
  {
    string filename = path + "T.txt";
    readIntFromASCIIFile(filename, T);
  }
  {
    string filename = path + "N_RWG.txt";
    readIntFromASCIIFile(filename, E);
  }
  {
    string filename = path + "S.txt";
    readIntFromASCIIFile(filename, S);
  }
  // we resize the arrays
  vertexes_coord.resize(V, 3);
  triangles_surfaces.resize(T);
  isClosedSurface.resize(S);
  cubes_RWGsNumbers.resize(E);
  cubes_centroids.resize(C, 3);
  RWGNumber_signedTriangles.resize(E, 2);
  RWGNumber_edgeVertexes.resize(E, 2);
  RWGNumber_oppVertexes.resize(E, 2);
  RWGNumber_CFIE_OK.resize(E);
  // reading RWG_numbers
  {
    string filename = path + "vertexes_coord.txt";
    readDoubleBlitzArray2DFromBinaryFile(filename, vertexes_coord);
  }
  {
    string filename = path + "triangles_surfaces.txt";
    readIntBlitzArray1DFromBinaryFile(filename, triangles_surfaces);
  }
  {
    string filename = path + "isClosedSurface.txt";
    readIntBlitzArray1DFromBinaryFile(filename, isClosedSurface);
  }
  {
    string filename = path + "cubes_RWGsNumbers.txt";
    readIntBlitzArray1DFromBinaryFile(filename, cubes_RWGsNumbers);
  }
  {
    string filename = path + "cubes_centroids.txt";
    readDoubleBlitzArray2DFromBinaryFile(filename, cubes_centroids);
  }
  {
    string filename = path + "RWGNumber_signedTriangles.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_signedTriangles);
  }
  {
    string filename = path + "RWGNumber_edgeVertexes.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_edgeVertexes);
  }
  {
    string filename = path + "RWGNumber_oppVertexes.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_oppVertexes);
  }
  {
    string filename = path + "RWGNumber_CFIE_OK.txt";
    readIntBlitzArray1DFromBinaryFile(filename, RWGNumber_CFIE_OK);
  }

};

void Mesh::resizeToZero() {
  triangles_surfaces.resize(0);
  isClosedSurface.resize(0);
  vertexes_coord.resize(0,0);
  cubes_RWGsNumbers.resize(0);
  cubes_centroids.resize(0,0);
  RWGNumber_signedTriangles.resize(0,0);
  RWGNumber_edgeVertexes.resize(0,0);
  RWGNumber_oppVertexes.resize(0,0);
  RWGNumber_CFIE_OK.resize(0);
}
*/
/******************************************************/
/******************** end of Mesh *********************/
/******************************************************/

/******************************************************/
/********************** LocalMesh *********************/
/******************************************************/

LocalMesh::LocalMesh(void){};

LocalMesh::LocalMesh(const string path) {setLocalMeshFromFile(path);};

void LocalMesh::copyLocalMesh(const LocalMesh& localMeshToCopy) /// copy member function
{
  N_local_RWG = localMeshToCopy.N_local_RWG;
  localRWGNumbers.resize(localMeshToCopy.localRWGNumbers.size());
  localRWGNumbers = localMeshToCopy.localRWGNumbers;
  reallyLocalRWGNumbers.resize(localMeshToCopy.reallyLocalRWGNumbers.size());
  reallyLocalRWGNumbers = localMeshToCopy.reallyLocalRWGNumbers;
  localRWGNumber_CFIE_OK.resize(localMeshToCopy.localRWGNumber_CFIE_OK.size());
  localRWGNumber_CFIE_OK = localMeshToCopy.localRWGNumber_CFIE_OK;
  localRWGNumber_signedTriangles.resize(localMeshToCopy.localRWGNumber_signedTriangles.extent(0), localMeshToCopy.localRWGNumber_signedTriangles.extent(1));
  localRWGNumber_signedTriangles = localMeshToCopy.localRWGNumber_signedTriangles;
  localRWGNumber_trianglesCoord.resize(localMeshToCopy.localRWGNumber_trianglesCoord.extent(0), localMeshToCopy.localRWGNumber_trianglesCoord.extent(1));
  localRWGNumber_trianglesCoord = localMeshToCopy.localRWGNumber_trianglesCoord;
}

LocalMesh::LocalMesh(const LocalMesh& localMeshToCopy) /// copy constructor
{
  copyLocalMesh(localMeshToCopy);
}

LocalMesh& LocalMesh::operator=(const LocalMesh& localMeshToCopy) { /// copy assignment
  copyLocalMesh(localMeshToCopy);
  return *this;
}

LocalMesh::~LocalMesh() {
  resizeToZero();
  localRWGNumbers.free();
  reallyLocalRWGNumbers.free();
  localRWGNumber_CFIE_OK.free();
  localRWGNumber_signedTriangles.free();
  localRWGNumber_trianglesCoord.free();
}

void LocalMesh::resizeToZero() {
  localRWGNumbers.resize(0);
  reallyLocalRWGNumbers.resize(0);
  localRWGNumber_CFIE_OK.resize(0);
  localRWGNumber_signedTriangles.resize(0,0);
  localRWGNumber_trianglesCoord.resize(0,0);
}

void LocalMesh::setLocalMeshFromFile(const string path)
{
  {
    string filename = path + "N_local_RWG.txt";
    readIntFromASCIIFile(filename, N_local_RWG);
  }
  // we resize the arrays
  localRWGNumbers.resize(N_local_RWG);
  reallyLocalRWGNumbers.resize(N_local_RWG);
  localRWGNumber_CFIE_OK.resize(N_local_RWG);
  localRWGNumber_signedTriangles.resize(N_local_RWG, 2);
  localRWGNumber_trianglesCoord.resize(N_local_RWG, 12);
  // reading RWG_numbers
  {
    string filename = path + "localRWGNumbers.txt";
    readIntBlitzArray1DFromBinaryFile(filename, localRWGNumbers);
  }
  for (int i=0 ; i<N_local_RWG ; ++i) reallyLocalRWGNumbers(i) = i;
  {
    string filename = path + "localRWGNumber_CFIE_OK.txt";
    readIntBlitzArray1DFromBinaryFile(filename, localRWGNumber_CFIE_OK);
  }
  {
    string filename = path + "localRWGNumber_signedTriangles.txt";
    readIntBlitzArray2DFromBinaryFile(filename, localRWGNumber_signedTriangles);
  }
  {
    string filename = path + "localRWGNumber_trianglesCoord.txt";
    readFloatBlitzArray2DFromBinaryFile(filename, localRWGNumber_trianglesCoord);
  }
};

void LocalMesh::writeLocalMeshToFile(const string path)
{
  {
    string filename = path + "N_local_RWG.txt";
    writeIntToASCIIFile(filename, N_local_RWG);
  }
  // writing the arrays
  {
    string filename = path + "localRWGNumbers.txt";
    writeIntBlitzArray1DToBinaryFile(filename, localRWGNumbers);
  }
  {
    string filename = path + "localRWGNumber_CFIE_OK.txt";
    writeIntBlitzArray1DToBinaryFile(filename, localRWGNumber_CFIE_OK);
  }
  {
    string filename = path + "localRWGNumber_signedTriangles.txt";
    writeIntBlitzArray2DToBinaryFile(filename, localRWGNumber_signedTriangles);
  }
  {
    string filename = path + "localRWGNumber_trianglesCoord.txt";
    writeFloatBlitzArray2DToBinaryFile(filename, localRWGNumber_trianglesCoord);
  }
};

/******************************************************/
/***************** end of LocalMesh *******************/
/******************************************************/

/******************************************************/
/************** various mesh functions ****************/
/******************************************************/

void compute_indexesEqualEdges(std::vector<std::vector<int> >& indexesEqualEdges,
                               const blitz::Array<int, 1>& indexesEqualPreceding,
                               const blitz::Array<int, 1>& ind_sorted_e_v)
{
  indexesEqualEdges.reserve(indexesEqualPreceding.size());
  int index = 0;
  for (unsigned int i=0 ; i<indexesEqualPreceding.size() ; i++) {
    const int indexEqualPreceding = indexesEqualPreceding(i);
    std::vector<int> tmp;
    tmp.resize(2);
    tmp[0] = ind_sorted_e_v(indexEqualPreceding-1);
    tmp[1] = ind_sorted_e_v(indexEqualPreceding);
    if (index>0) {
      if (tmp[0]==indexesEqualEdges[index-1].back()) indexesEqualEdges[index-1].push_back(tmp[1]);
      else {
        indexesEqualEdges.push_back(tmp);
        index += 1;
      }
    }
    else {
      indexesEqualEdges.push_back(tmp);
      index += 1;
    }
  }
  for (unsigned int i=0 ; i<indexesEqualEdges.size() ; i++) std::vector<int> (indexesEqualEdges[i]).swap(indexesEqualEdges[i]);
}

void compute_edgeNumber_vertexes(blitz::Array<int, 2>& edgeNumber_vertexes,
                                 const std::vector<std::vector<int> >& indexesEqualEdges,
                                 const blitz::Array<int, 2>& col_sorted_e_v)
{
  // renumbering of the edges
  int N_equalEdges = indexesEqualEdges.size();
  edgeNumber_vertexes.resize(N_equalEdges, 2);
  for (int i=0 ; i<N_equalEdges ; i++) {
    for (int j=0 ; j<2; j++) edgeNumber_vertexes(i, j) = col_sorted_e_v(indexesEqualEdges[i][0], j);
  }
}

void compute_edgeNumber_triangles(blitz::Array<int, 2>& edgeNumber_triangles,
                                  const std::vector<std::vector<int> >& indexesEqualEdges)
{
  unsigned int N_columns = indexesEqualEdges[0].size();
  for (unsigned int i=1 ; i<indexesEqualEdges.size() ; i++) {
    if (indexesEqualEdges[i].size() > N_columns) N_columns = indexesEqualEdges[i].size();
  }
  edgeNumber_triangles.resize(indexesEqualEdges.size(), N_columns);
  edgeNumber_triangles = -1;
  for (unsigned int i=0 ; i<indexesEqualEdges.size() ; i++) {
    for (unsigned int j=0 ; j<indexesEqualEdges[i].size() ; j++) edgeNumber_triangles(int(i), int(j)) = indexesEqualEdges[i][j]/3;
  }
}

void compute_triangle_adjacentTriangles(std::vector<std::vector<int> >& triangle_adjacentTriangles,
                                        const blitz::Array<int, 2>& edgeNumber_triangles,
                                        const int T)
{
  triangle_adjacentTriangles.resize(T);
  const int N_edges = edgeNumber_triangles.extent(0), N_columns_edges = edgeNumber_triangles.extent(1);
  for (int i=0 ; i<N_edges ; i++) {
    std::vector<int> adjacent_triangles;
    for (int j=0 ; j<N_columns_edges ; j++) {
      if (edgeNumber_triangles(i, j) == -1) break;
      else adjacent_triangles.push_back(edgeNumber_triangles(i, j));
    }
    const int N_adj_triangles = adjacent_triangles.size();
    const bool IS_JUNCTION = (N_adj_triangles>2);
    for (int j=0 ; j<N_adj_triangles ; j++) {
      const int tn = adjacent_triangles[j];
      for (int k=0 ; k<N_adj_triangles ; k++) {
        const int tk = adjacent_triangles[k];
        if ( (tk!=tn) && (IS_JUNCTION) ) triangle_adjacentTriangles[tn].push_back(-tk-1);
        else if (tk!=tn) triangle_adjacentTriangles[tn].push_back(tk);
      }
    }
  }
  // memory enhancement...
  for (int i=0 ; i<T ; i++) std::vector<int> (triangle_adjacentTriangles[i]).swap(triangle_adjacentTriangles[i]);
}

/*****************************/
/* reorder_triangle_vertexes */
/*****************************/

void change_triangle_circulation(blitz::Array<int, 2>& triangle_vertexes, const int t0, const int t1)
/**
 *  This function reorders the nodes of triangle t1 such that the common edge
 *  between reference triangle t0 and t1 is parcouru following one direction and
 *  then following the opposite. In this way we will have compatible normals
 */
{
  // nodes of reference triangle
  const int node00 = triangle_vertexes(t0, 0);
  const int node01 = triangle_vertexes(t0, 1);
  const int node02 = triangle_vertexes(t0, 2);
  // nodes of triangle to change
  const int node10 = triangle_vertexes(t1, 0);
  const int node11 = triangle_vertexes(t1, 1);
  const int node12 = triangle_vertexes(t1, 2);
  // circulation on these triangles is a succession of paths between nodes: 0 -> 1, 1 -> 2, 2 -> 0
  blitz::Array<int, 2> t0_circ(3, 2), t1_circ(3, 2);
  t0_circ = node00, node01,
            node01, node02,
            node02, node00;
  t1_circ = node10, node11,
            node11, node12,
            node12, node10;
  // "coincide_circ": a variable that will tell if the edge is parcouru alike by each triangle.
  // In this case, one will have to reorder the nodes of t1.
  bool coincide_circ = false;
  for (int i=0 ; i<3 ; i++) {
    if (coincide_circ) break;
    for (int j=0 ; j<3 ; j++) {
      if ( (t0_circ(i, 0)==t1_circ(j, 0)) && (t0_circ(i, 1)==t1_circ(j, 1)) ) {
        coincide_circ = true;
        break;
      }
    }
  }
  if (coincide_circ) {
    triangle_vertexes(t1, 1) = node12;
    triangle_vertexes(t1, 2) = node11;
  }
}

bool array_has_zero(int & index, const blitz::Array<int, 1>& x)
/* this function checks if an int array has zero and where */
{
  const int N = x.size();
  index = 0;
  bool HAS_ZERO = false;
  for (index=0 ; index<N ; index++) {
    if (x(index)==0) {
      HAS_ZERO = true;
      break;
    }
  }
  return HAS_ZERO;
}

void compute_list_t_to_reorder(std::vector<int>& list_t_to_reorder,
                               std::vector<int>& list_calling_t_to_reorder,
                               const int calling_t,
                               const std::vector<std::vector<int> >& triangle_adjacentTriangles,
                               const blitz::Array<int, 1>& is_triangle_reordered,
                               blitz::Array<int, 1>& is_triangle_in_list)
{
  std::vector<int> calling_t_adjacentTriangles = triangle_adjacentTriangles[calling_t];
  // we fill in list_t_to_reorder
  for (unsigned int i=0 ; i<calling_t_adjacentTriangles.size() ; i++) {
    int tn = calling_t_adjacentTriangles[i];
    bool is_triangle_not_adjacent_via_junction = true;
    if (tn < 0) { // then the two triangles are adjacent via junction
      tn = -tn - 1;
      is_triangle_not_adjacent_via_junction = false;
    }
    if ( (is_triangle_reordered(tn)==0) && (is_triangle_in_list(tn)==0) && (is_triangle_not_adjacent_via_junction) ) {
      is_triangle_in_list(tn) = 1;
      list_t_to_reorder.push_back(tn);     
    }
  }
  // creation of the list_calling_t
  list_calling_t_to_reorder.resize(list_t_to_reorder.size());
  for (unsigned int i=0 ; i<list_calling_t_to_reorder.size() ; i++) list_calling_t_to_reorder[i] = calling_t;
}

void reorder_triangle_vertexes(blitz::Array<int, 2>& triangle_vertexes,
                               blitz::Array<int, 1>& triangles_surfaces,
                               const blitz::Array<double, 2>& vertexes_coord,
                               const std::vector<std::vector<int> >& triangle_adjacentTriangles)
{
  const int T = triangle_adjacentTriangles.size();
  // triangles_surfaces initialization
  triangles_surfaces.resize(T);
  triangles_surfaces = -1;
  // is_triangle_reordered tells if a triangle has been reordered (1) or not (0)
  blitz::Array<int, 1> is_triangle_reordered(T), is_triangle_in_list(T);
  std::vector<int> is_triangle_in_listTmp;
  is_triangle_reordered = 0;
  is_triangle_in_list = 0;
  int surface_number = -1, index_first_zero = 0;
  while (array_has_zero(index_first_zero, is_triangle_reordered)) {
    int t_start = index_first_zero;
    surface_number++;
    is_triangle_reordered(t_start) = 1;
    triangles_surfaces(t_start) = surface_number;
    // construction of the list of triangles that will be reordered
    std::vector<int> list_t_to_reorder, list_calling_t;
    compute_list_t_to_reorder(list_t_to_reorder, list_calling_t, t_start, triangle_adjacentTriangles, is_triangle_reordered, is_triangle_in_list);
    while (list_t_to_reorder.size()>0) {
      int t = list_t_to_reorder.back();
      int calling_t = list_calling_t.back();
      is_triangle_in_listTmp.push_back(t);
      list_t_to_reorder.pop_back();
      list_calling_t.pop_back();
      change_triangle_circulation(triangle_vertexes, calling_t, t);
      is_triangle_reordered(t) = 1;
      is_triangle_reordered(calling_t) = 1;
      triangles_surfaces(t) = surface_number;
      triangles_surfaces(calling_t) = surface_number;
      // we now augment list_t_to_reorder and list_calling_t
      std::vector<int> list_t_to_reorderTmp, list_calling_tTmp;
      compute_list_t_to_reorder(list_t_to_reorderTmp, list_calling_tTmp, t, triangle_adjacentTriangles, is_triangle_reordered, is_triangle_in_list);
      for (unsigned int i=0 ; i<list_t_to_reorderTmp.size() ; i++) {
        list_t_to_reorder.push_back(list_t_to_reorderTmp[i]);
        list_calling_t.push_back(list_calling_tTmp[i]);
      }
    } // end inner while
  } // end outer while

  // We now reorient the triangles so that the normals point outward for each surface.
  // For this, we find the "highest" triangle of the surface, i.e. the one that has
  // the highest centroid and we see if its normal points "upward". If it is not the
  // case, we simply swap the middle and last columns of "triangle_vertexes" for triangles on s.
  const int S = surface_number;
  blitz::Array<double, 1> triangles_centroids_z(T);
  for (int i=0 ; i<T ; i++) {
    double centroid_z = 0.0;
    for (int j=0 ; j<3 ; j++) centroid_z += vertexes_coord(triangle_vertexes(i, j), 2);
    triangles_centroids_z(i) = centroid_z/3.0;
  }

  for (int s=0 ; s<S+1 ; s++) {
    std::vector<int> ind_t_on_s;
    ind_t_on_s.reserve(T);
    // let's get the indexes of triangles on surface s
    for (int j=0 ; j<T ; j++) {
      if (triangles_surfaces(j)==s) ind_t_on_s.push_back(j);
    }
    // now we want the index of the one that has highest centroid
    double max_height_centroids_s = triangles_centroids_z(ind_t_on_s[0]);
    int index_max_height_centroids_s = ind_t_on_s[0];
    for (unsigned int j=1 ; j<ind_t_on_s.size() ; j++) {
      if (triangles_centroids_z(ind_t_on_s[j])>max_height_centroids_s) {
        max_height_centroids_s = triangles_centroids_z(ind_t_on_s[j]);
        index_max_height_centroids_s = ind_t_on_s[j];
      }
    }
    // let's now see if the normal points "upward"
    blitz::Range all = blitz::Range::all();
    const blitz::Array<double, 1> r0(vertexes_coord(triangle_vertexes(index_max_height_centroids_s, 0), all));
    const blitz::Array<double, 1> r1(vertexes_coord(triangle_vertexes(index_max_height_centroids_s, 1), all));
    const blitz::Array<double, 1> r2(vertexes_coord(triangle_vertexes(index_max_height_centroids_s, 2), all));
    const blitz::Array<double, 1> r1_r0(r1-r0);
    const blitz::Array<double, 1> r2_r0(r2-r0);
    const double triangle_normal_z = r1_r0(0) * r2_r0(1) - r1_r0(1) * r2_r0(0);
    // if (triangle_normal_z<0) we have to swap the columns of triangle_vertexes
    if (triangle_normal_z<0) {
      for (unsigned int j=0 ; j<ind_t_on_s.size() ; j++) {
        const int index = ind_t_on_s[j];
        const int v1 = triangle_vertexes(index, 1), v2 = triangle_vertexes(index, 2);
        triangle_vertexes(index, 1) = v2;
        triangle_vertexes(index, 2) = v1;
      }
    }
  }
}

/*****************************/
/*    is_surface_closed      */
/*****************************/

void is_surface_closed(blitz::Array<int, 1>& is_closed_surface,
                       blitz::Array<std::vector<int>, 2>& connected_surfaces,
                       blitz::Array<std::vector<int>, 2>& potential_closed_surfaces,
                       const blitz::Array<int, 1>& triangles_surfaces,
                       const blitz::Array<int, 2>& edgeNumber_triangles)
/**
 *  This function determines if a surface (made up of triangles) is closed
 *  or not. It also provides relationships between surfaces (linked or not).
 */
{
  const int S = max(triangles_surfaces) + 1, T = triangles_surfaces.size();
  is_closed_surface.resize(S);
  connected_surfaces.resize(S, S);
  potential_closed_surfaces.resize(S, S);
  // we count the number of trianges on each surface
  blitz::Array<int, 1> NUMBER_TRIANGLES_IN_SURFACE(S);
  NUMBER_TRIANGLES_IN_SURFACE = 0;
  for (int j=0 ; j<T ; j++) NUMBER_TRIANGLES_IN_SURFACE(triangles_surfaces(j)) += 1;
  // we now count the number of INNER edges for each surface.
  // the edges that are junctions receive a special treatment:
  // only if the edge has two triangles on the given surface,
  // can it be counted as an inner edge, which will then be
  // counted in NUMBER_EDGES_IN_SURFACE
  blitz::Array<int, 1> NUMBER_EDGES_IN_SURFACE(S);
  NUMBER_EDGES_IN_SURFACE = 0;
  const int N_edges = edgeNumber_triangles.extent(0), N_columns_edges = edgeNumber_triangles.extent(1);
  for (int i=0 ; i<N_edges ; i++) {
    blitz::Array<int, 1> surfaces_appeared_already(S);
    surfaces_appeared_already = 0;
    std::vector<int> triangles_tmp;
    for (int j=0 ; j<N_columns_edges ; j++) {
      if (edgeNumber_triangles(i, j) == -1) break;
      else triangles_tmp.push_back(edgeNumber_triangles(i, j));
    }
    // if (triangles_tmp.size()>2) we have a junction
    if (triangles_tmp.size()>2) {
      for (unsigned int j=0 ; j<triangles_tmp.size() ; j++) {
        const int t = triangles_tmp[j];
        const int surface = triangles_surfaces(t);
        if (surfaces_appeared_already(surface)==0) surfaces_appeared_already(surface) = 1;
        else NUMBER_EDGES_IN_SURFACE(surface) += 1;
      }
      std::vector<int> surfaces_present;
      for (unsigned int j=0 ; j<surfaces_appeared_already.size() ; j++) {
        if (surfaces_appeared_already(j)>0) surfaces_present.push_back(j);
      }
      // filling of connected_surfaces
      for (unsigned int index1=0 ; index1<surfaces_present.size() ; index1++) {
        for (unsigned int index2=index1+1 ; index2<surfaces_present.size() ; index2++) {
          const int s1 = min(surfaces_present[index1], surfaces_present[index2]);
          const int s2 = max(surfaces_present[index1], surfaces_present[index2]);
          connected_surfaces(s1, s2).push_back(i);
        }
      }

    } // end if (triangles_tmp.size()>2)
    else { // it is not a junction then
      const int surface = triangles_surfaces(triangles_tmp[0]);
      NUMBER_EDGES_IN_SURFACE(surface) += 1;
    }
  } // end for
  // filling of is_closed_surface
  for (int i=0 ; i<S ; i++) is_closed_surface(i) = ( (NUMBER_EDGES_IN_SURFACE(i) * 2)==(NUMBER_TRIANGLES_IN_SURFACE(i) * 3) ) * 1;
  // we now check for potential closed surfaces: combination of surfaces
  // which can be closed and on which we can therefore apply the CFIE
  for (int s0=0 ; s0<S ; s0++) {
    for (int s1=0 ; s1<S ; s1++) {
      const int N_connecting_edges = connected_surfaces(s0, s1).size();
      if (N_connecting_edges>0) {
        const int numberEdges0 = NUMBER_EDGES_IN_SURFACE(s0), numberEdges1 = NUMBER_EDGES_IN_SURFACE(s1);
        const int numberTriangles0 = NUMBER_TRIANGLES_IN_SURFACE(s0), numberTriangles1 = NUMBER_TRIANGLES_IN_SURFACE(s1);
        if ( ((numberEdges0 + numberEdges1 + N_connecting_edges)*2) == (3*(numberTriangles0 + numberTriangles1)) ) potential_closed_surfaces(s0, s1) = connected_surfaces(s0, s1);
      }
    }
  }
}



/*****************************************/
/* RWGNumber_signedTriangles_computation */
/*****************************************/

void RWGNumber_signedTriangles_computation(blitz::Array<int, 2>& RWGNumber_signedTriangles,
                                           blitz::Array<int, 2>& RWGNumber_edgeVertexes,
                                           const blitz::Array<int, 2>& edgeNumber_triangles,
                                           const blitz::Array<int, 2>& edgeNumber_vertexes,
                                           const blitz::Array<int, 1>& triangles_surfaces,
                                           const blitz::Array<int, 1>& is_closed_surface,
                                           const blitz::Array<int, 2>& triangle_vertexes,
                                           const blitz::Array<double, 2>& vertexes_coord)
{
  blitz::Range all = blitz::Range::all();
  const int N_edges = edgeNumber_triangles.extent(0), N_columns_edges = edgeNumber_triangles.extent(1);
  // RWGNumber_signedTrianglesTmp_1 is an array of fixed size N_edges. It will be equal to
  // edgeNumber_triangles if there are no junctions. If there are
  // junctions, RWGNumber_signedTrianglesTmp_2 will be non-empty. This is for limiting
  // the memory requirements on this part of the code.
  blitz::Array<int, 2> RWGNumber_signedTrianglesTmp(N_edges, 2);
  std::vector<std::vector<int> > RWGNumber_signedTrianglesTmp_1;
  std::vector<int> RWGNumber_edgeNumber; // needed for RWGNumber_edgeVertexes computation
  RWGNumber_edgeNumber.reserve(N_edges);
  for (int i=0 ; i<N_edges ; i++) RWGNumber_edgeNumber.push_back(i);
  // we now construct the RWGs. It will be done differently if we have a junction or not
  for (int i=0 ; i<N_edges ; i++) {
    std::vector<int> triangles;
    for (int j=0 ; j<N_columns_edges ; j++) {
      if (edgeNumber_triangles(i, j) == -1) break;
      else triangles.push_back(edgeNumber_triangles(i, j));
    }
    const int numberOfTriangles = triangles.size();
    // if (numberOfTriangles==2) we have a classic RWG
    if (numberOfTriangles==2) {
      RWGNumber_signedTrianglesTmp(i, 0) = triangles[0];
      RWGNumber_signedTrianglesTmp(i, 1) = triangles[1];
    }
    else { // else we have a junction
      // we compute the edge unit vector: zHatEdge
      const int n0 = edgeNumber_vertexes(i, 0), n1 = edgeNumber_vertexes(i, 1);
      blitz::Array<double, 1> r0(3), r1(3), r2(3);
      r0 = vertexes_coord(n0, all);
      r1 = vertexes_coord(n1, all);
      const blitz::Array<double, 1> r1_r0(r1 - r0);
      const blitz::Array<double, 1> zHatEdge(r1_r0 / sqrt(sum(r1_r0 * r1_r0)));
      blitz::Array<double, 1> xHatEdge(3), yHatEdge(3);
      blitz::Array<double, 1> triangles_angles(numberOfTriangles);
      // construction of the edge local coordinate system: it is based on the first triangle
      // this is for sorting the triangles wrt their respective angles
      // because we cannot have a RWG which is bisected by another RWG
      for (int j=0 ; j<numberOfTriangles ; j++) {
        const int tr = triangles[j];
        const blitz::Array<int, 1> tr_nodes(triangle_vertexes(tr, all));
        // for each triangle we find the corner opposite to the edge
        for (unsigned int n2=0 ; n2<tr_nodes.size() ; n2++) {
          if ( (tr_nodes(n2)!=n0) && (tr_nodes(n2)!=n1) ) {
            r2 = vertexes_coord(tr_nodes(n2), all);
          }
        }
        const blitz::Array<double, 1> r2_r0(r2 - r0);
        const blitz::Array<double, 1> r2_r0_perpendicular(r2_r0 - sum(r2_r0 * zHatEdge) * zHatEdge);
        const blitz::Array<double, 1> r2_r0_perpendicularHat(r2_r0_perpendicular/sqrt(sum(r2_r0_perpendicular * r2_r0_perpendicular)));
        // construction of the edge local coordinate system
        if (j==0) {// the first triangle will be the reference: will define the xHat vector
          xHatEdge = r2_r0_perpendicularHat;
          // yHat = zHat X xHat
          yHatEdge(0) = zHatEdge(1) * xHatEdge(2) - zHatEdge(2) * xHatEdge(1);
          yHatEdge(1) = zHatEdge(2) * xHatEdge(0) - zHatEdge(0) * xHatEdge(2);
          yHatEdge(2) = zHatEdge(0) * xHatEdge(1) - zHatEdge(1) * xHatEdge(0);
          triangles_angles(0) = 0.0;
        }
        else {
          const double Cos = sum(r2_r0_perpendicularHat * xHatEdge);
          const double Sin = sum(r2_r0_perpendicularHat * yHatEdge);
          if (Sin>=0) triangles_angles(j) = acos(Cos);
          else triangles_angles(j) = 2.0*M_PI - acos(Cos);
        }
      }
      // we now sort the triangles wrt their respective position wrt xHatEdge
      // we need an argsort type function, given by the Dictionary class (see mesh.h)
      std::vector< Dictionary<double, int> > triangles_anglesToIndexes;
      triangles_anglesToIndexes.reserve(numberOfTriangles);
      for (int j=0 ; j<numberOfTriangles ; j++) triangles_anglesToIndexes.push_back(Dictionary<double, int> (triangles_angles(j), j));
      sort(triangles_anglesToIndexes.begin(), triangles_anglesToIndexes.end());
      std::vector<int> sortedTriangles, sortedTrianglesSurfaces;
      sortedTriangles.reserve(numberOfTriangles);
      sortedTrianglesSurfaces.reserve(numberOfTriangles);
      for (int j=0 ; j<numberOfTriangles ; j++) {
        sortedTriangles.push_back(triangles[triangles_anglesToIndexes[j].getVal()]);
        sortedTrianglesSurfaces.push_back(triangles_surfaces(sortedTriangles.back()));
      }
      // we now form all the possible RWGs with the sorted triangles.
      // Normally none of these RWGs can be bisected by a triangle now
      // THE FOLLOWING ONLY WORKS FOR METAL-METAL JUNCTIONS
      std::vector<std::vector<int> >possibleTrianglesPairsForRWGs;
      for (int j=1 ; j<numberOfTriangles ; j++) {
        if (sortedTrianglesSurfaces[j] != sortedTrianglesSurfaces[j-1]) {
          std::vector<int> temp;
          temp.push_back(sortedTriangles[j-1]);
          temp.push_back(sortedTriangles[j]);
          possibleTrianglesPairsForRWGs.push_back(temp);
        }
      }
      // the last possibility of RWG: between the last and first triangle
      if (possibleTrianglesPairsForRWGs.size()<numberOfTriangles-1) { // if our list is too small
        if (sortedTrianglesSurfaces[numberOfTriangles-1] != sortedTrianglesSurfaces[0]) {
          std::vector<int> temp;
          temp.push_back(sortedTriangles[numberOfTriangles-1]);
          temp.push_back(sortedTriangles[0]);
          possibleTrianglesPairsForRWGs.push_back(temp);
        }
      }

      for (unsigned int j=1 ; j<possibleTrianglesPairsForRWGs.size() ; j++) {
        RWGNumber_signedTrianglesTmp_1.push_back(possibleTrianglesPairsForRWGs[j]);
        RWGNumber_edgeNumber.push_back(i);
      }
      // current line of RWGNumber_signedTrianglesTmp_1
      RWGNumber_signedTrianglesTmp(i, 0) = possibleTrianglesPairsForRWGs[0][0];
      RWGNumber_signedTrianglesTmp(i, 1) = possibleTrianglesPairsForRWGs[0][1];
    }
  }
  // we can finally construct the final RWGNumber_signedTriangles array
  const int N_RWG = N_edges + RWGNumber_signedTrianglesTmp_1.size();
  RWGNumber_signedTriangles.resize(N_RWG, 2);
  for (int i=0; i<N_RWG ; i++) {
    if (i<N_edges) {
      RWGNumber_signedTriangles(i, 0) = RWGNumber_signedTrianglesTmp(i, 0);
      RWGNumber_signedTriangles(i, 1) = RWGNumber_signedTrianglesTmp(i, 1);
    }
    else {
      RWGNumber_signedTriangles(i, 0) = RWGNumber_signedTrianglesTmp_1[i-N_edges][0];
      RWGNumber_signedTriangles(i, 1) = RWGNumber_signedTrianglesTmp_1[i-N_edges][1];
    }
  }
  RWGNumber_signedTrianglesTmp_1.clear();
  RWGNumber_signedTrianglesTmp.free();
  RWGNumber_edgeVertexes.resize(N_RWG, 2);
  for (int i=0 ; i<N_RWG ; i++) {
    // computation of RWGNumber_edgeVertexes
    const int t0 = RWGNumber_signedTriangles(i, 0), t1 = RWGNumber_signedTriangles(i, 1);
    const int s0 = triangles_surfaces(t0), s1 = triangles_surfaces(t1);
    int t = t0;
    if (s0 != s1) {
      // we should choose the triangle which belongs to a closed surface as the reference triangle
      if (is_closed_surface(s0)==1) t = t0;
      else if (is_closed_surface(s1)==1) t = t1;
      else t = t0;
    }
    const int e0 = edgeNumber_vertexes(RWGNumber_edgeNumber[i], 0), e1 = edgeNumber_vertexes(RWGNumber_edgeNumber[i], 1);
    const int n0 = triangle_vertexes(t, 0), n1 = triangle_vertexes(t, 1), n2 = triangle_vertexes(t, 2);
    if ( ((e0==n0) && (e1==n1)) || ((e0==n1) && (e1==n2)) || ((e0==n2) && (e1==n0)) ) {
      RWGNumber_edgeVertexes(i, 0) = e0;
      RWGNumber_edgeVertexes(i, 1) = e1;
    }
    else {
      RWGNumber_edgeVertexes(i, 0) = e1;
      RWGNumber_edgeVertexes(i, 1) = e0;
    }
  }
}


void RWGNumber_oppVertexes_computation(blitz::Array<int, 2>& RWGNumber_oppVertexes,
                                       const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                                       const blitz::Array<int, 2>& RWGNumber_edgeVertexes,
                                       const blitz::Array<int, 2>& triangle_vertexes)
{
  const int N_RWG = RWGNumber_signedTriangles.extent(0);
  RWGNumber_oppVertexes.resize(N_RWG, 2);
  for (int i=0 ; i<N_RWG ; i++) {
    const int e0 = RWGNumber_edgeVertexes(i, 0), e1 = RWGNumber_edgeVertexes(i, 1);
    for (int j=0 ; j<2 ; j++) {
      const int t = RWGNumber_signedTriangles(i, j);
      const int node0 = triangle_vertexes(t, 0),
                node1 = triangle_vertexes(t, 1),
                node2 = triangle_vertexes(t, 2);
      if ( (node0!=e0) && (node0!=e1) ) RWGNumber_oppVertexes(i, j) = node0;
      else if ( (node1!=e0) && (node1!=e1) ) RWGNumber_oppVertexes(i, j) = node1;
      else if ( (node2!=e0) && (node2!=e1) ) RWGNumber_oppVertexes(i, j) = node2;
      else {
        cout << "mesh.cpp::RWGNumber_oppVertexes_computation: error in opp nodes computation, RWGNumber = " << i << ", j = " << j << endl;
        exit(1);
      }
    }
  }
}

