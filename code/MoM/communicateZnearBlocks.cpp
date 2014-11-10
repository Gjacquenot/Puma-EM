#include <fstream>
#include <iostream>
#include <string>
#include <blitz/array.h>
#include <mpi.h>

using namespace blitz;

#include "GetMemUsage.h"
#include "readWriteBlitzArrayFromFile.h"

/****************************************************************************/
/******************************* main ***************************************/
/****************************************************************************/

int main(int argc, char* argv[]) {

  MPI::Init();
  int my_id = MPI::COMM_WORLD.Get_rank();
  int num_procs = MPI::COMM_WORLD.Get_size();

  string simuDir = ".";
  if ( argc > 2 ) {
     if( string(argv[1]) == "--simudir" ) simuDir = argv[2];
  }

  MPI_Status irecv_status, isend_status_1, irecv_status_1;
  MPI_Request isend_request_1, irecv_request_1;
  
  //  general variables
  const string Z_BLOCKS_PATH = simuDir + "/tmp" + intToString(my_id) + "/Z_tmp/";
  int itemsize;
  readIntFromASCIIFile(Z_BLOCKS_PATH + "itemsize.txt", itemsize);

  if (my_id==0) cout << "Exchanging Z_near blocks for preconditioner construction......" << endl;
  flush(cout);
  MPI_Barrier(MPI_COMM_WORLD);

  // we loop on the process numbers
  for (int recv_proc_number = 0; recv_proc_number<num_procs; recv_proc_number++) {
    blitz::Array<int, 1> CubesNumbersToReceive, ChunkNumbersToReceive;

    for (int send_proc_number = 0; send_proc_number<num_procs; send_proc_number++) {
      if (recv_proc_number != send_proc_number) {
        int N_cubes;
        blitz::Array<int, 1> CubesNumbersToSend, ChunkNumbersToSend;
        if (my_id==send_proc_number) {
          readIntBlitzArray1DFromASCIIFile(Z_BLOCKS_PATH + "CubesNumbersToSendToP" + intToString(recv_proc_number) + ".txt", CubesNumbersToSend);
          readIntBlitzArray1DFromASCIIFile(Z_BLOCKS_PATH + "ChunkNumbersToSendToP" + intToString(recv_proc_number) + ".txt", ChunkNumbersToSend);
          N_cubes = CubesNumbersToSend.size();
          MPI_Send(&N_cubes, 1, MPI_INT, recv_proc_number, recv_proc_number, MPI_COMM_WORLD);
          MPI_Send(CubesNumbersToSend.data(), CubesNumbersToSend.size(), MPI_INT, recv_proc_number, recv_proc_number+1, MPI_COMM_WORLD);
          MPI_Send(ChunkNumbersToSend.data(), ChunkNumbersToSend.size(), MPI_INT, recv_proc_number, recv_proc_number+2, MPI_COMM_WORLD);
        }
        if (my_id==recv_proc_number) {
          MPI_Recv(&N_cubes, 1, MPI_INT, send_proc_number, my_id, MPI_COMM_WORLD, &irecv_status);
          CubesNumbersToReceive.resize(N_cubes);
          ChunkNumbersToReceive.resize(N_cubes);
          MPI_Recv(CubesNumbersToReceive.data(), CubesNumbersToReceive.size(), MPI_INT, send_proc_number, my_id+1, MPI_COMM_WORLD, &irecv_status);
          MPI_Recv(ChunkNumbersToReceive.data(), ChunkNumbersToReceive.size(), MPI_INT, send_proc_number, my_id+2, MPI_COMM_WORLD, &irecv_status);
        }
        const int NCubesToReceive = CubesNumbersToReceive.size(), NCubesToSend = CubesNumbersToSend.size();
        blitz::Array< MPI_Request, 1 > Array_isend_request(NCubesToSend), Array_irecv_request(NCubesToReceive);
        blitz::Array< MPI_Status, 1 > Array_isend_status(NCubesToSend), Array_irecv_status(NCubesToReceive);

        // we first send the cubes mesh data
        // 1) we first have to read the int arrays
        blitz::Array<int, 1> N_IntArraysToReceive(NCubesToReceive), N_IntArraysToSend(NCubesToSend);
        blitz::Array<blitz::Array<int, 1>, 1> IntArraysToReceive(NCubesToReceive), IntArraysToSend(NCubesToSend);
        if (my_id == send_proc_number)  {
          for (int i=0 ; i<NCubesToSend ; ++i) {
            const string CHUNK_PATH = "chunk" + intToString(ChunkNumbersToSend(i)), 
                         filename_N_Int = intToString(CubesNumbersToSend(i)) + "_N_IntArrays.txt";
            const string File_N_IntArray = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename_N_Int;
            // now reading the arrays themselves
            const string filename_IntArray = intToString(CubesNumbersToSend(i)) + "_IntArrays.txt";
            const string File_IntArray = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename_IntArray;
            // reading the size of the int array
            blitz::ifstream ifs(File_IntArray.c_str(), blitz::ios::binary);
            ifs.seekg (0, blitz::ios::end);
            int length = ifs.tellg();
            ifs.close();
            N_IntArraysToSend(i) = length/4;

            IntArraysToSend(i).resize(N_IntArraysToSend(i));
            readIntBlitzArray1DFromBinaryFile(File_IntArray, IntArraysToSend(i));
          }
          // we now send the sizes arrays to the receiving process
          MPI_Isend(N_IntArraysToSend.data(), N_IntArraysToSend.size(), MPI_INT, recv_proc_number, 23, MPI_COMM_WORLD, &isend_request_1);
        }
        // the receives of the arrays sizes
        if (my_id==recv_proc_number) {
          MPI_Irecv(N_IntArraysToReceive.data(), N_IntArraysToReceive.size(), MPI_INT, send_proc_number, 23, MPI_COMM_WORLD, &irecv_request_1);
        }
        if (my_id == recv_proc_number) MPI_Wait(&irecv_request_1, &irecv_status_1);
        if (my_id == send_proc_number) MPI_Wait(&isend_request_1, &isend_status_1);

        // we can now post the arrays themselves
        if (my_id == send_proc_number)  {
          for (int i=0 ; i<NCubesToSend ; ++i) {
             int cubeNumber = CubesNumbersToSend(i);
             MPI_Isend(IntArraysToSend(i).data(), IntArraysToSend(i).size(), MPI_INT, recv_proc_number, cubeNumber, MPI_COMM_WORLD, &Array_isend_request(i));
          }
        }
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int cubeNumber = CubesNumbersToReceive(i);
            IntArraysToReceive(i).resize(N_IntArraysToReceive(i));
            MPI_Irecv(IntArraysToReceive(i).data(), IntArraysToReceive(i).size(), MPI_INT, send_proc_number, cubeNumber, MPI_COMM_WORLD, &Array_irecv_request(i));
          }
        }
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) MPI_Wait(&Array_irecv_request(i), &Array_irecv_status(i));
        }
        if (my_id == send_proc_number) {
          for (int i=0 ; i<NCubesToSend ; ++i) MPI_Wait(&Array_isend_request(i), &Array_isend_status(i));
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // finally we write the communicated files to disk
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int src = send_proc_number;
            if (src == send_proc_number) {
              const string CHUNK_PATH = "chunk" + intToString(ChunkNumbersToReceive(i)); 
              string fileToWrite = Z_BLOCKS_PATH + CHUNK_PATH + "/" + intToString(CubesNumbersToReceive(i)) + "_IntArrays.txt";
              writeIntBlitzArray1DToBinaryFile(fileToWrite, IntArraysToReceive(i));
            }
          }
        }

        // we then send the Z matrices
        blitz::Array<blitz::Array<std::complex<float>, 2>, 1> send_buff(NCubesToSend);
        blitz::Array<blitz::Array<std::complex<float>, 2>, 1> recv_buff(NCubesToReceive);
        if (my_id == send_proc_number)  {
          for (int i=0 ; i<NCubesToSend ; ++i) {
            int cubeNumber = CubesNumbersToSend(i), chunkNumber = ChunkNumbersToSend(i), dest = recv_proc_number;
            if (my_id == dest) cout << "PROCESS " << my_id <<  ": SENDING FAILED!!!!!!!" << endl;
            const int Nl = IntArraysToSend(i)(0), Nc = IntArraysToSend(i)(1);
            send_buff(i).resize(Nl, Nc);
            const string CHUNK_PATH = "chunk" + intToString(chunkNumber), filename = intToString(cubeNumber);
            const string fileToRead = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename;
            if (itemsize==8) readComplexFloatBlitzArray2DFromBinaryFile(fileToRead, send_buff(i));
            int BUF_SIZE = send_buff(i).size();
            MPI_Isend(send_buff(i).data(), BUF_SIZE, MPI_COMPLEX, dest, cubeNumber, MPI_COMM_WORLD, &Array_isend_request(i));
          }
        }
        // we post the receives of recv_proc_number
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int cubeNumber = CubesNumbersToReceive(i), src = send_proc_number;
            if (my_id == src) cout << "RECEIVING FAILED!!!!!!!" << endl;
            const int Nl = IntArraysToReceive(i)(0), Nc = IntArraysToReceive(i)(1);
            recv_buff(i).resize(Nl, Nc);
            // we post the receives only if my_id == recv_proc_number, in order to not post
            // all the receives at the same time, which could deadlock
            // because of too many communications
            int BUF_SIZE = recv_buff(i).size();
            MPI_Irecv(recv_buff(i).data(), BUF_SIZE, MPI_COMPLEX, src, cubeNumber, MPI_COMM_WORLD, &Array_irecv_request(i));
          }
        }

        //if (my_id == recv_proc_number) cout << "Receiving Process " << recv_proc_number << ", Sending Process " << send_proc_number << ": All commands posted. Waiting to achieve......" << endl;
        // we wait for all the communications to finish
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) MPI_Wait(&Array_irecv_request(i), &Array_irecv_status(i));
        }
        if (my_id == send_proc_number) {
          for (int i=0 ; i<NCubesToSend ; ++i) MPI_Wait(&Array_isend_request(i), &Array_isend_status(i));
        }
        //if (my_id == recv_proc_number) cout << "Receiving Process " << recv_proc_number << ", Sending Process " << send_proc_number << ": Achieved!!!......" << endl;

        // finally we write the communicated files to disk
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int cubeNumber = CubesNumbersToReceive(i), chunkNumber = ChunkNumbersToReceive(i), src = send_proc_number;
            if (src == send_proc_number) {
              const string CHUNK_PATH = "chunk" + intToString(chunkNumber), filename = intToString(cubeNumber);
              const string fileToWrite = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename;
              blitz::ofstream fout(fileToWrite.c_str(), blitz::ios::binary);
              fout.write((char *)(recv_buff(i).data()), recv_buff(i).size()*itemsize);
              fout.close();
            }
          }
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;

  MPI::Finalize();
  return 0;
}

