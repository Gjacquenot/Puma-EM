#include "scatter_mesh_per_cube.h"

void list_testSrc_RWGs_computation(blitz::Array<int, 1>& test_RWGs_numbers,
                                   blitz::Array<int, 1>& src_RWGs_numbers,
                                   blitz::Array<int, 1>& cube_neighborsIndexes,
                                   const int cubeNumber,
                                   const blitz::Array<int, 1>& cubes_NeighborsIndexes,
                                   const blitz::Array<int, 1>& cubes_N_neighbors, 
                                   const blitz::Array<int, 1>& cube_startIndex_neighbors,
                                   const blitz::Array<int, 1>& cubes_RWGsNumbers, 
                                   const blitz::Array<int, 1>& cube_N_RWGs, 
                                   const blitz::Array<int, 1>& cube_startIndex_RWGs)
{
  // first the test RWGs
  {
    const int N_test = cube_N_RWGs(cubeNumber);
    test_RWGs_numbers.resize(N_test);
    int startIndex = cube_startIndex_RWGs(cubeNumber);
    for (int i=0; i<N_test; i++) test_RWGs_numbers(i) = cubes_RWGsNumbers(i + startIndex);
  }
  
  // then the src RWGs
  // we construct an array containing the neighbors indexes
  const int N_neighbors = cubes_N_neighbors(cubeNumber);
  cube_neighborsIndexes.resize(N_neighbors);
  {
    int startIndex = cube_startIndex_neighbors(cubeNumber);
    for (int i=0; i<N_neighbors; i++) cube_neighborsIndexes(i) = cubes_NeighborsIndexes(i + startIndex);
  }
  // we compute the size of src_RWGs_numbers
  int N_src = 0;
  for (int i=0; i<N_neighbors; i++) {
    const int neighborNumber = cube_neighborsIndexes(i);
    N_src += cube_N_RWGs(neighborNumber);
  }
  src_RWGs_numbers.resize(N_src);
  // now we fill in src_RWGs_numbers
  int index = 0;
  for (int i=0; i<N_neighbors; i++) {
    const int neighborNumber = cube_neighborsIndexes(i);
    const int N_RWGs = cube_N_RWGs(neighborNumber);
    const int startIndex = cube_startIndex_RWGs(neighborNumber);
    for (int j=0; j<N_RWGs; j++) {
      src_RWGs_numbers(index) = cubes_RWGsNumbers(j + startIndex);
      index++;
    }
  }
}

void compute_cube_arrays_from_mesh(blitz::Array<int, 1>& cubeIntArrays,
                                   blitz::Array<double, 1>& cubeDoubleArrays,
                                   const int cubeNumber,
                                   const int S,
                                   const blitz::Array<int, 1>& cubes_NeighborsIndexes, // sizeInt approx C*27
                                   const blitz::Array<int, 1>& cubes_N_neighbors, // sizeInt C
                                   const blitz::Array<int, 1>& cube_startIndex_neighbors, // sizeInt C
                                   const blitz::Array<int, 1>& cubes_RWGsNumbers, // sizeInt N_RWG
                                   const blitz::Array<int, 1>& cube_N_RWGs, // sizeInt C
                                   const blitz::Array<int, 1>& cube_startIndex_RWGs, // sizeInt C
                                   const blitz::Array<int, 1>& RWGNumber_CFIE_OK, // sizeInt N_RWG
                                   const blitz::Array<int, 1>& RWGNumber_M_CURRENT_OK, // sizeInt N_RWG
                                   const blitz::Array<int, 2>& RWGNumber_signedTriangles, // sizeInt (N_RWG, 2)
                                   const blitz::Array<int, 2>& RWGNumber_edgeVertexes, // sizeInt (N_RWG, 2)
                                   const blitz::Array<int, 2>& RWGNumber_oppVertexes, // sizeInt (N_RWG, 2)
                                   const blitz::Array<double, 2>& vertexes_coord, // sizeDouble (V, 3)
                                   const blitz::Array<double, 2>& cubes_centroids) // sizeDouble (C, 3)
{
  blitz::Range all = blitz::Range::all();
  blitz::Array<int, 1> test_RWGs_numbers, src_RWGs_numbers, cube_neighborsIndexes;
  list_testSrc_RWGs_computation(test_RWGs_numbers, src_RWGs_numbers, cube_neighborsIndexes, cubeNumber, cubes_NeighborsIndexes, cubes_N_neighbors, cube_startIndex_neighbors, cubes_RWGsNumbers, cube_N_RWGs, cube_startIndex_RWGs);
  const int N_RWG_test = test_RWGs_numbers.size();
  const int N_RWG_src = src_RWGs_numbers.size();
  const int N_neighbors = cube_neighborsIndexes.size();

  cubeIntArrays.resize(5 + N_RWG_src + 2*N_RWG_src + 4*N_RWG_src + N_RWG_test + N_RWG_src + N_neighbors);
  blitz::Array<int, 1> localTestSrcRWGNumber_signedTriangles(N_RWG_src * 2);
  for (int i=0; i<N_RWG_src; i++) {
    int index = src_RWGs_numbers(i);
    localTestSrcRWGNumber_signedTriangles(2*i) = RWGNumber_signedTriangles(index, 0);
    localTestSrcRWGNumber_signedTriangles(2*i + 1) = RWGNumber_signedTriangles(index, 1);
  }

  blitz::Array<int, 1> localSrcRWGNumber_M_CURRENT_OK(N_RWG_src);
  for (int i=0; i<N_RWG_src; i++) localSrcRWGNumber_M_CURRENT_OK(i) = RWGNumber_M_CURRENT_OK(src_RWGs_numbers(i));

  blitz::Array<int, 1> localTestRWGNumber_CFIE_OK(N_RWG_test);
  for (int i=0; i<N_RWG_test; i++) localTestRWGNumber_CFIE_OK(i) = RWGNumber_CFIE_OK(test_RWGs_numbers(i));

  blitz::Array<int, 1> localTestSrcRWGNumber_nodes(N_RWG_src * 4);
  for (int i=0; i<N_RWG_src; i++) {
    const int RWG_number = src_RWGs_numbers(i);
    const int nodeEdge_0 = RWGNumber_edgeVertexes(RWG_number, 0);
    const int nodeEdge_1 = RWGNumber_edgeVertexes(RWG_number, 1);
    const int oppNodeEdge_0 = RWGNumber_oppVertexes(RWG_number, 0);
    const int oppNodeEdge_1 = RWGNumber_oppVertexes(RWG_number, 1);
    localTestSrcRWGNumber_nodes(4*i) = oppNodeEdge_0;
    localTestSrcRWGNumber_nodes(4*i + 1) = nodeEdge_0;
    localTestSrcRWGNumber_nodes(4*i + 2) = nodeEdge_1;
    localTestSrcRWGNumber_nodes(4*i + 3) = oppNodeEdge_1;
  }
  // we create a dictionary that we will need for creating the local_vertexes_coord
  std::vector< Dictionary<int, int> > oldNodeNumber_to_index;
  oldNodeNumber_to_index.reserve(N_RWG_src * 4);
  for (int i=0; i<N_RWG_src * 4; i++) oldNodeNumber_to_index.push_back( Dictionary<int, int> (localTestSrcRWGNumber_nodes(i), i) );

  // we sort the dictionary by its "localTestSrcRWGNumber_nodes" values (the "key" field)
  sort(oldNodeNumber_to_index.begin(), oldNodeNumber_to_index.end());

  // we now count the number of different nodes
  int N_nodes = 1;
  for (int i=1; i<N_RWG_src * 4; i++) {
    if (oldNodeNumber_to_index[i].getKey() != oldNodeNumber_to_index[i-1].getKey()) N_nodes += 1;
  }

  blitz::Array<double, 2> local_vertexes_coord(N_nodes + 1, 3); // the "+1" is for the cube centroid
  // we now start replacing the nodes numbers in localTestSrcRWGNumber_nodes by their new, local number
  int newNodeNumber = 0;
  int index = oldNodeNumber_to_index[0].getVal();
  localTestSrcRWGNumber_nodes(index) = newNodeNumber;
  local_vertexes_coord(newNodeNumber, all) = vertexes_coord(oldNodeNumber_to_index[0].getKey(), all);
  for (int i=1; i<N_RWG_src * 4; i++) {
    if (oldNodeNumber_to_index[i].getKey() != oldNodeNumber_to_index[i-1].getKey()) {
      newNodeNumber++;
      local_vertexes_coord(newNodeNumber, all) = vertexes_coord(oldNodeNumber_to_index[i].getKey(), all);
    }
    index = oldNodeNumber_to_index[i].getVal();
    localTestSrcRWGNumber_nodes(index) = newNodeNumber;
  }
  local_vertexes_coord(N_nodes, all) = cubes_centroids(cubeNumber, all);

  // now we fill up the cubeIntArrays
  int startIndex = 5;
  int stopIndex = startIndex + N_RWG_src;
  for (int i=startIndex; i<stopIndex; i++) cubeIntArrays(i) = src_RWGs_numbers(i-startIndex);

  startIndex = stopIndex;
  stopIndex = startIndex + 2*N_RWG_src;
  for (int i=startIndex; i<stopIndex; i++) cubeIntArrays(i) = localTestSrcRWGNumber_signedTriangles(i-startIndex);

  startIndex = stopIndex;
  stopIndex = startIndex + 4*N_RWG_src;
  for (int i=startIndex; i<stopIndex; i++) cubeIntArrays(i) = localTestSrcRWGNumber_nodes(i-startIndex);

  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_test;
  for (int i=startIndex; i<stopIndex; i++) cubeIntArrays(i) = localTestRWGNumber_CFIE_OK(i-startIndex);

  startIndex = stopIndex;
  stopIndex = startIndex + N_RWG_src;
  for (int i=startIndex; i<stopIndex; i++) cubeIntArrays(i) = localSrcRWGNumber_M_CURRENT_OK(i-startIndex);

  startIndex = stopIndex;
  stopIndex = startIndex + N_neighbors;
  for (int i=startIndex; i<stopIndex; i++) cubeIntArrays(i) = cube_neighborsIndexes(i-startIndex);

  // final moves
  cubeIntArrays(0) = N_RWG_test;
  cubeIntArrays(1) = N_RWG_src;
  cubeIntArrays(2) = N_neighbors;
  cubeIntArrays(3) = N_nodes;
  cubeIntArrays(4) = S;

  cubeDoubleArrays.resize( (N_nodes + 1) * 3 );
  for (int i=0; i<N_nodes+1; i++) {
    cubeDoubleArrays(3*i) = local_vertexes_coord(i, 0);
    cubeDoubleArrays(3*i+1) = local_vertexes_coord(i, 1);
    cubeDoubleArrays(3*i+2) = local_vertexes_coord(i, 2);    
  }
}

/*void Mg_listsOfZnearBlocks_ToTransmitAndReceive(const blitz::Array<int, 1>& local_chunkNumber_to_cubesNumbers,
                                                const blitz::Array<int, 1>& local_cubeNumber_to_chunkNumbers,
                                                const blitz::Array<blitz::Array<int, 1>, 1>& allCubeIntArrays,
                                                const string Z_TMP_DATA_PATH)
{
  const int N_local_cubes = local_chunkNumber_to_cubesNumbers.size();
  const int num_procs = MPI::COMM_WORLD.Get_size();
  const int my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0;
  int ierror;

  // we will need a multimap
  multimap <int, int> local_cubes, neighbor_cubes;
  multimap <int, int>::iterator itr, itr_2;
  for (int i=0; i<N_local_cubes; i++) local_cubes.insert(pair<int, int> (local_chunkNumber_to_cubesNumbers(i), local_cubeNumber_to_chunkNumbers(i)));

  // we construct a list of the non local cubes neighbors, that is,
  // the neighbors which are contained in another process
  std::vector<int> listOfNonLocalCubesNeighborsTmp;
  for (int i=0; i<N_local_cubes; i++) {
    const int N_neighbors = allCubeIntArrays(i)(2);
    const int N_int = allCubeIntArrays(i).size();
    blitz::Array<int, 1> neighbors(N_neighbors);
    const int startIndex = N_int-N_neighbors;
    for (int j=0; j<N_neighbors; j++) neighbors(j) = allCubeIntArrays(i)(startIndex + j);
    for (int j=0; j<N_neighbors; j++) {
      itr = neighbor_cubes.find(neighbors(j));
      if ( itr == neighbor_cubes.end() ) {
        neighbor_cubes.insert( pair<int, int> (neighbors(j), 0) );
        itr_2 = local_cubes.find(neighbors(j));
        if ( itr_2 == local_cubes.end() ) listOfNonLocalCubesNeighborsTmp.push_back(neighbors(j));
      }
    }
  }
  sort(listOfNonLocalCubesNeighborsTmp.begin(), listOfNonLocalCubesNeighborsTmp.end());
  // we put the non local cubes into an array for MPI sending
  const int N_non_local = listOfNonLocalCubesNeighborsTmp.size();
  blitz::Array<int, 1> listOfNonLocalCubesNeighbors(N_non_local);
  for (int i=0; i<N_non_local; i++) listOfNonLocalCubesNeighbors(i) = listOfNonLocalCubesNeighborsTmp[i];

  for (int P=0; P<num_procs; P++) {
    // first we broadcast the listOfNonLocalCubesNeighbors from process P to all processes
    int N_non_local_cubes_process_P;
    if (my_id==P) N_non_local_cubes_process_P = N_non_local;
    MPI_Bcast(&N_non_local_cubes_process_P, 1, MPI_INT, P, MPI_COMM_WORLD);
    blitz::Array<int, 1> listOfNonLocalCubesNeighbors_process_P(N_non_local_cubes_process_P);
    if (my_id==P) {
      for (int i=0; i<N_non_local_cubes_process_P; i++) listOfNonLocalCubesNeighbors_process_P(i) = listOfNonLocalCubesNeighbors(i);
    }
    MPI_Bcast(listOfNonLocalCubesNeighbors_process_P.data(), N_non_local_cubes_process_P, MPI_INT, P, MPI_COMM_WORLD);
    
    // now we evaluate for each process the list of cubes it has to send to process P
    blitz::Array<int, 1> listCubesToSendToProcessP(0), listChunksToSendToProcessP(0);
    int N_to_send = 0;
    if (my_id!=P) {
      std::vector<int> listCubesToSendToProcessP_tmp, listChunksToSendToProcessP_tmp;
      for (int i=0; i<N_non_local_cubes_process_P; i++) {
        const int Non_local_cube_number = listOfNonLocalCubesNeighbors_process_P(i);
        itr = local_cubes.find(Non_local_cube_number);
        if ( itr != local_cubes.end() ) {
          listCubesToSendToProcessP_tmp.push_back(Non_local_cube_number);
          listChunksToSendToProcessP_tmp.push_back((*itr).second);
        }
      }
      N_to_send = listCubesToSendToProcessP_tmp.size();
      listCubesToSendToProcessP.resize(N_to_send);
      listChunksToSendToProcessP.resize(N_to_send);
      for (int i=0; i<N_to_send; i++) {
        listCubesToSendToProcessP(i) = listCubesToSendToProcessP_tmp[i];
        listChunksToSendToProcessP(i) = listChunksToSendToProcessP_tmp[i];
      }
      // now we write the arrays to the disk
      string filename = Z_TMP_DATA_PATH + "CubesNumbersToSendToP" + intToString(P) + ".txt";
      writeIntBlitzArray1DToASCIIFile(filename, listCubesToSendToProcessP);
      filename = Z_TMP_DATA_PATH + "ChunkNumbersToSendToP" + intToString(P) + ".txt";
      writeIntBlitzArray1DToASCIIFile(filename, listChunksToSendToProcessP);
    } 
  }
}*/

void scatter_mesh_per_cube(blitz::Array<blitz::Array<int, 1>, 1>& allCubeIntArrays,
                           blitz::Array<blitz::Array<double, 1>, 1>& allCubeDoubleArrays,
                           const string simuDir,
                           const blitz::Array<int, 1>& local_ChunksNumbers,
                           const blitz::Array<int, 1>& local_chunkNumber_N_cubesNumbers,
                           const blitz::Array<int, 1>& local_chunkNumber_to_cubesNumbers,
                           const blitz::Array<int, 1>& local_cubeNumber_to_chunkNumbers)
{
  int ierror;
  const int num_procs = MPI::COMM_WORLD.Get_size();
  const int my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0;
  MPI_Status status;

  // general variables
  const string SIMU_DIR = simuDir;
  const string TMP = SIMU_DIR + "/tmp" + intToString(my_id);
  const string MESH_DATA_PATH = TMP + "/mesh/";
  const string Z_TMP_DATA_PATH = TMP + "/Z_tmp/";
  string filename;

  // reading mesh data
  int C, N_RWG, S, V;
  blitz::Array<int, 1> cubes_RWGsNumbers, cube_N_RWGs, cube_startIndex_RWGs;
  blitz::Array<int, 1> cubes_neighborsIndexes, cube_N_neighbors, cube_startIndex_neighbors;
  blitz::Array<int, 1> RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK;
  blitz::Array<int, 2> RWGNumber_edgeVertexes, RWGNumber_oppVertexes, RWGNumber_signedTriangles;
  blitz::Array<double, 2> cubes_centroids, vertexes_coord;

  if (my_id==0)
  {
    // reading the arrays sizes
    filename = MESH_DATA_PATH + "C.txt";
    readIntFromASCIIFile(filename, C);
    
    filename = MESH_DATA_PATH + "N_RWG.txt";
    readIntFromASCIIFile(filename, N_RWG);
    
    filename = MESH_DATA_PATH + "S.txt";
    readIntFromASCIIFile(filename, S);

    filename = MESH_DATA_PATH + "V.txt";
    readIntFromASCIIFile(filename, V);

    // resizing
    vertexes_coord.resize(V, 3);
    cube_N_RWGs.resize(C);
    cube_startIndex_RWGs.resize(C);
    cube_N_neighbors.resize(C);
    cube_startIndex_neighbors.resize(C);
    cubes_RWGsNumbers.resize(N_RWG);
    cubes_centroids.resize(C, 3);
    RWGNumber_edgeVertexes.resize(N_RWG, 2);
    RWGNumber_oppVertexes.resize(N_RWG, 2);
    RWGNumber_signedTriangles.resize(N_RWG, 2);
    RWGNumber_CFIE_OK.resize(N_RWG);
    RWGNumber_M_CURRENT_OK.resize(N_RWG);

    // reading the arrays
    filename = MESH_DATA_PATH + "cubes_centroids.txt";
    readDoubleBlitzArray2DFromBinaryFile(filename, cubes_centroids);
    
    filename = MESH_DATA_PATH + "cubes_RWGsNumbers.txt";
    readIntBlitzArray1DFromBinaryFile(filename, cubes_RWGsNumbers);
    
    filename = MESH_DATA_PATH + "cube_N_RWGs.txt";
    readIntBlitzArray1DFromBinaryFile(filename, cube_N_RWGs);
    
    // creation of cube_startIndex_RWGs
    int startIndex = 0;
    for (int i=0; i<C; i++) {
      cube_startIndex_RWGs(i) = startIndex;
      startIndex += cube_N_RWGs(i);
    }
    
    filename = MESH_DATA_PATH + "cube_N_neighbors.txt";
    readIntBlitzArray1DFromBinaryFile(filename, cube_N_neighbors);

    // creation of cube_startIndex_neighbors
    startIndex = 0;
    for (int i=0; i<C; i++) {
      cube_startIndex_neighbors(i) = startIndex;
      startIndex += cube_N_neighbors(i);
    }
    int N_neighbors = startIndex;
    // startIndex now represents the total number of neighbors
    cubes_neighborsIndexes.resize(N_neighbors);

    filename = MESH_DATA_PATH + "cubes_neighborsIndexes.txt";
    readIntBlitzArray1DFromBinaryFile(filename, cubes_neighborsIndexes);

    filename = MESH_DATA_PATH + "vertexes_coord.txt";
    readDoubleBlitzArray2DFromBinaryFile(filename, vertexes_coord);

    filename = MESH_DATA_PATH + "RWGNumber_edgeVertexes.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_edgeVertexes);

    filename = MESH_DATA_PATH + "RWGNumber_oppVertexes.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_oppVertexes);

    filename = MESH_DATA_PATH + "RWGNumber_signedTriangles.txt";
    readIntBlitzArray2DFromBinaryFile(filename, RWGNumber_signedTriangles);

    filename = MESH_DATA_PATH + "RWGNumber_CFIE_OK.txt";
    readIntBlitzArray1DFromBinaryFile(filename, RWGNumber_CFIE_OK);

    filename = MESH_DATA_PATH + "RWGNumber_M_CURRENT_OK.txt";
    readIntBlitzArray1DFromBinaryFile(filename, RWGNumber_M_CURRENT_OK);
  }
  MPI_Bcast(&N_RWG, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&C, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&S, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // now we have to read the local chunks numbers, and the local cubes numbers. 
  // They will be communicated to the master process. Copy beginning of mpi_mlfma.cpp for this.
  int N_local_Chunks = local_ChunksNumbers.size();
  int N_local_cubes = local_chunkNumber_to_cubesNumbers.size();

  // preparing the terrain for gathering the cubes numbers for each process  
  blitz::Array<int, 1> NumberOfCubesPerProcess;
  if (my_id==0) NumberOfCubesPerProcess.resize(num_procs);
  ierror = MPI_Gather(&N_local_cubes, 1, MPI_INT, NumberOfCubesPerProcess.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  // we prepare the scounts and displacements for gathering the cubes numbers
  blitz::Array<int, 1> MPI_Gatherv_scounts, MPI_Gatherv_displs; // it only matters for process 0
  if (my_id==0) {
    MPI_Gatherv_scounts.resize(num_procs);
    MPI_Gatherv_displs.resize(num_procs);
    int displacement = 0;
    for (int i=0 ; i<num_procs ; ++i) {
      MPI_Gatherv_scounts(i) = NumberOfCubesPerProcess(i);
      MPI_Gatherv_displs(i) = displacement;
      displacement += MPI_Gatherv_scounts(i);
    }
  }

  blitz::Array<int, 1> process_cubesNumbers;
  if (my_id==0) process_cubesNumbers.resize(C);
  blitz::Array<int, 1> local_cubes(local_chunkNumber_to_cubesNumbers);
  ierror = MPI_Gatherv(local_cubes.data(), N_local_cubes, MPI_INT, process_cubesNumbers.data(), MPI_Gatherv_scounts.data(), MPI_Gatherv_displs.data(), MPI_INT, 0,  MPI_COMM_WORLD);
  
  blitz::Array<int, 1> cubeArraysSizes(2);
  allCubeIntArrays.resize(N_local_cubes);
  allCubeDoubleArrays.resize(N_local_cubes);
  if (my_id==master) {
    for (int receive_id=num_procs-1; receive_id>-1; receive_id--) {
      // first we find the RWGs for each cube of process receive_id
      const int startIndexOfCube = MPI_Gatherv_displs(receive_id);
      const int Ncubes = MPI_Gatherv_scounts(receive_id);
      for (int i=0; i<Ncubes; i++) {
        blitz::Array<int, 1> cubeIntArrays;
        blitz::Array<double, 1> cubeDoubleArrays; 
        const int cubeNumber = process_cubesNumbers(startIndexOfCube + i);
        compute_cube_arrays_from_mesh(cubeIntArrays, cubeDoubleArrays, cubeNumber, S, cubes_neighborsIndexes, cube_N_neighbors, cube_startIndex_neighbors, cubes_RWGsNumbers, cube_N_RWGs, cube_startIndex_RWGs, RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK, RWGNumber_signedTriangles, RWGNumber_edgeVertexes, RWGNumber_oppVertexes, vertexes_coord, cubes_centroids);
        cubeArraysSizes(0) = cubeIntArrays.size();
        cubeArraysSizes(1) = cubeDoubleArrays.size();
        if (receive_id!=master) {
          // we send the arrays to the receiving process
          MPI_Send(cubeArraysSizes.data(), cubeArraysSizes.size(), MPI_INT, receive_id, receive_id, MPI_COMM_WORLD);
          MPI_Send(cubeIntArrays.data(), cubeIntArrays.size(), MPI_INT, receive_id, cubeNumber, MPI_COMM_WORLD);
          MPI_Send(cubeDoubleArrays.data(), cubeDoubleArrays.size(), MPI_DOUBLE, receive_id, cubeNumber+1, MPI_COMM_WORLD);
        }
        else {
          // we are on master node, we write the cube arrays to disk
          allCubeIntArrays(i).resize(cubeArraysSizes(0));
          allCubeDoubleArrays(i).resize(cubeArraysSizes(1));
          allCubeIntArrays(i) = cubeIntArrays;
          allCubeDoubleArrays(i) = cubeDoubleArrays;
        }
      }
    }
  }
  if (my_id!=master) {
    // we receive the cubes arrays and write them to disk
    for (int i=0; i<N_local_cubes; i++) {
      const int cubeNumber = local_chunkNumber_to_cubesNumbers(i);
      MPI_Recv(cubeArraysSizes.data(), cubeArraysSizes.size(), MPI_INT, 0, my_id, MPI_COMM_WORLD, &status);
      allCubeIntArrays(i).resize(cubeArraysSizes(0));
      allCubeDoubleArrays(i).resize(cubeArraysSizes(1));
      MPI_Recv(allCubeIntArrays(i).data(), allCubeIntArrays(i).size(), MPI_INT, 0, cubeNumber, MPI_COMM_WORLD, &status);
      MPI_Recv(allCubeDoubleArrays(i).data(), allCubeDoubleArrays(i).size(), MPI_DOUBLE, 0, cubeNumber+1, MPI_COMM_WORLD, &status);
    }
  }
}

