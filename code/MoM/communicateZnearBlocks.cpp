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

  MPI_Status isend_status, irecv_status, isend_status_1, irecv_status_1;
  MPI_Request isend_request, irecv_request, isend_request_1, irecv_request_1;
  
  //  general variables
  const string Z_BLOCKS_PATH = simuDir + "/tmp" + intToString(my_id) + "/Z_tmp/";
  int itemsize;
  readIntFromASCIIFile(Z_BLOCKS_PATH + "itemsize.txt", itemsize);

  if (my_id==0) cout << "Exchanging Z_near blocks for preconditioner construction......" << endl;
  flush(cout);
  int ierror = MPI_Barrier(MPI::COMM_WORLD), ierror_1;

  // we loop on the process numbers
  for (int recv_proc_number = 0; recv_proc_number<num_procs; recv_proc_number++) {
    blitz::Array<int, 1> CubesNumbersToReceive, ChunkNumbersToReceive;
    blitz::Array<int, 2> dimensionsOfCubesToReceive;

    for (int send_proc_number = 0; send_proc_number<num_procs; send_proc_number++) {
      if (recv_proc_number != send_proc_number) {
        blitz::Array<int, 1> CubesNumbersToSend, ChunkNumbersToSend;
        blitz::Array<int, 2> dimensionsOfCubesToSend;
        if (my_id==recv_proc_number) {
          readIntBlitzArray1DFromASCIIFile(Z_BLOCKS_PATH + "CubesNumbersToReceiveFromP" + intToString(send_proc_number) + ".txt", CubesNumbersToReceive);
          readIntBlitzArray1DFromASCIIFile(Z_BLOCKS_PATH + "ChunkNumbersToReceiveFromP" + intToString(send_proc_number) + ".txt", ChunkNumbersToReceive);
          dimensionsOfCubesToReceive.resize(CubesNumbersToReceive.size(), 2);
          ierror = MPI_Irecv(dimensionsOfCubesToReceive.data(), dimensionsOfCubesToReceive.size(), MPI::INT, send_proc_number, 22, MPI::COMM_WORLD, &irecv_request);
        }
        if (my_id==send_proc_number) {
          readIntBlitzArray1DFromASCIIFile(Z_BLOCKS_PATH + "CubesNumbersToSendToP" + intToString(recv_proc_number) + ".txt", CubesNumbersToSend);
          readIntBlitzArray1DFromASCIIFile(Z_BLOCKS_PATH + "ChunkNumbersToSendToP" + intToString(recv_proc_number) + ".txt", ChunkNumbersToSend);
          readIntBlitzArray2DFromASCIIFile(Z_BLOCKS_PATH + "dimensionsOfCubesToSendToP" + intToString(recv_proc_number) + ".txt", dimensionsOfCubesToSend);
          ierror = MPI_Isend(dimensionsOfCubesToSend.data(), dimensionsOfCubesToSend.size(), MPI::INT, recv_proc_number, 22, MPI::COMM_WORLD, &isend_request);
        }
        if (my_id == recv_proc_number) ierror = MPI_Wait(&irecv_request, &irecv_status);
        if (my_id == send_proc_number) ierror = MPI_Wait(&isend_request, &isend_status);

        const int NCubesToReceive = CubesNumbersToReceive.size(), NCubesToSend = CubesNumbersToSend.size();
        blitz::Array< MPI_Request, 1 > Array_isend_request(NCubesToSend), Array_irecv_request(NCubesToReceive);
        blitz::Array< MPI_Status, 1 > Array_isend_status(NCubesToSend), Array_irecv_status(NCubesToReceive);

        // we first send the cubes mesh data
        // 1) we first have to read the int arrays and double arrays
        blitz::Array<int, 1> N_IntArraysToReceive(NCubesToReceive), N_IntArraysToSend(NCubesToSend);
        blitz::Array<int, 1> N_DoubleArraysToReceive(NCubesToReceive), N_DoubleArraysToSend(NCubesToSend);
        blitz::Array<blitz::Array<int, 1>, 1> IntArraysToReceive(NCubesToReceive), IntArraysToSend(NCubesToSend);
        blitz::Array<blitz::Array<double, 1>, 1> DoubleArraysToReceive(NCubesToReceive), DoubleArraysToSend(NCubesToSend);
        if (my_id == send_proc_number)  {
          for (int i=0 ; i<NCubesToSend ; ++i) {
            const string CHUNK_PATH = "chunk" + intToString(ChunkNumbersToSend(i)), 
                         filename_N_Int = intToString(CubesNumbersToSend(i)) + "_N_IntArrays.txt";
            const string File_N_IntArray = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename_N_Int;
            // now reading the arrays themselves
            const string filename_DoubleArray = intToString(CubesNumbersToSend(i)) + "_DoubleArrays.txt", 
                         filename_IntArray = intToString(CubesNumbersToSend(i)) + "_IntArrays.txt";
            const string File_DoubleArray = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename_DoubleArray, 
                         File_IntArray = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename_IntArray;
            // reading the size of the int array
            blitz::ifstream ifs(File_IntArray.c_str(), blitz::ios::binary);
            ifs.seekg (0, blitz::ios::end);
            int length = ifs.tellg();
            ifs.close();
            N_IntArraysToSend(i) = length/4;

            IntArraysToSend(i).resize(N_IntArraysToSend(i));
            readIntBlitzArray1DFromBinaryFile(File_IntArray, IntArraysToSend(i));
            // Now the arrays of Double
            N_DoubleArraysToSend(i) = IntArraysToSend(i)(3) * 3 + 3;
            DoubleArraysToSend(i).resize(N_DoubleArraysToSend(i));
            readDoubleBlitzArray1DFromBinaryFile(File_DoubleArray, DoubleArraysToSend(i));
          }
          // we now send the sizes arrays to the receiving process
          ierror = MPI_Isend(N_DoubleArraysToSend.data(), N_DoubleArraysToSend.size(), MPI::INT, recv_proc_number, 22, MPI::COMM_WORLD, &isend_request);
          ierror_1 = MPI_Isend(N_IntArraysToSend.data(), N_IntArraysToSend.size(), MPI::INT, recv_proc_number, 23, MPI::COMM_WORLD, &isend_request_1);
        }
        // the receives of the arrays sizes
        if (my_id==recv_proc_number) {
          ierror = MPI_Irecv(N_DoubleArraysToReceive.data(), N_DoubleArraysToReceive.size(), MPI::INT, send_proc_number, 22, MPI::COMM_WORLD, &irecv_request);
          ierror_1 = MPI_Irecv(N_IntArraysToReceive.data(), N_IntArraysToReceive.size(), MPI::INT, send_proc_number, 23, MPI::COMM_WORLD, &irecv_request_1);
        }
        if (my_id == recv_proc_number) ierror = MPI_Wait(&irecv_request, &irecv_status);
        if (my_id == recv_proc_number) ierror_1 = MPI_Wait(&irecv_request_1, &irecv_status_1);
        if (my_id == send_proc_number) ierror = MPI_Wait(&isend_request, &isend_status);
        if (my_id == send_proc_number) ierror_1 = MPI_Wait(&isend_request_1, &isend_status_1);

        // we can now post the arrays themselves
        if (my_id == send_proc_number)  {
          for (int i=0 ; i<NCubesToSend ; ++i) {
             int cubeNumber = CubesNumbersToSend(i);
             ierror = MPI_Isend(IntArraysToSend(i).data(), IntArraysToSend(i).size(), MPI::INT, recv_proc_number, cubeNumber, MPI::COMM_WORLD, &Array_isend_request(i));
          }
        }
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int cubeNumber = CubesNumbersToReceive(i);
            IntArraysToReceive(i).resize(N_IntArraysToReceive(i));
            ierror = MPI_Irecv(IntArraysToReceive(i).data(), IntArraysToReceive(i).size(), MPI::INT, send_proc_number, cubeNumber, MPI::COMM_WORLD, &Array_irecv_request(i));
          }
        }
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) ierror = MPI_Wait(&Array_irecv_request(i), &Array_irecv_status(i));
        }
        if (my_id == send_proc_number) {
          for (int i=0 ; i<NCubesToSend ; ++i) ierror = MPI_Wait(&Array_isend_request(i), &Array_isend_status(i));
        }
        ierror = MPI_Barrier(MPI::COMM_WORLD);

        // the doubles arrays
        if (my_id == send_proc_number)  {
          for (int i=0 ; i<NCubesToSend ; ++i) {
             int cubeNumber = CubesNumbersToSend(i);
             ierror = MPI_Isend(DoubleArraysToSend(i).data(), DoubleArraysToSend(i).size(), MPI::DOUBLE, recv_proc_number, cubeNumber, MPI::COMM_WORLD, &Array_isend_request(i));
          }
        }
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int cubeNumber = CubesNumbersToReceive(i);
            DoubleArraysToReceive(i).resize(N_DoubleArraysToReceive(i));
            ierror = MPI_Irecv(DoubleArraysToReceive(i).data(), DoubleArraysToReceive(i).size(), MPI::DOUBLE, send_proc_number, cubeNumber, MPI::COMM_WORLD, &Array_irecv_request(i));
          }
        }
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) ierror = MPI_Wait(&Array_irecv_request(i), &Array_irecv_status(i));
        }
        if (my_id == send_proc_number) {
          for (int i=0 ; i<NCubesToSend ; ++i) ierror = MPI_Wait(&Array_isend_request(i), &Array_isend_status(i));
        }
        ierror = MPI_Barrier(MPI::COMM_WORLD);

        // finally we write the communicated files to disk
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int src = send_proc_number;
            if (src == send_proc_number) {
              const string CHUNK_PATH = "chunk" + intToString(ChunkNumbersToReceive(i)); 
              {
                string fileToWrite = Z_BLOCKS_PATH + CHUNK_PATH + "/" + intToString(CubesNumbersToReceive(i)) + "_IntArrays.txt";
                writeIntBlitzArray1DToBinaryFile(fileToWrite, IntArraysToReceive(i));
              }
              {
                string fileToWrite = Z_BLOCKS_PATH + CHUNK_PATH + "/" + intToString(CubesNumbersToReceive(i)) + "_DoubleArrays.txt";
                writeDoubleBlitzArray1DToBinaryFile(fileToWrite, DoubleArraysToReceive(i));
              }

//              string fileToWrite = Z_BLOCKS_PATH + CHUNK_PATH + "/" + intToString(CubesNumbersToReceive(i)) + "_N_IntArrays.txt";
//              writeIntToASCIIFile(fileToWrite, N_IntArraysToReceive(i));
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
            const int Nl = dimensionsOfCubesToSend(i, 0), Nc = dimensionsOfCubesToSend(i, 1);
            send_buff(i).resize(Nl, Nc);
            const string CHUNK_PATH = "chunk" + intToString(chunkNumber), filename = intToString(cubeNumber);
            const string fileToRead = Z_BLOCKS_PATH + CHUNK_PATH + "/" + filename;
            if (itemsize==8) readComplexFloatBlitzArray2DFromBinaryFile(fileToRead, send_buff(i));
            int BUF_SIZE = send_buff(i).size();
            ierror = MPI_Isend(send_buff(i).data(), BUF_SIZE, MPI::COMPLEX, dest, cubeNumber, MPI::COMM_WORLD, &Array_isend_request(i));
          }
        }
        // we post the receives of recv_proc_number
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) {
            int cubeNumber = CubesNumbersToReceive(i), src = send_proc_number;
            if (my_id == src) cout << "RECEIVING FAILED!!!!!!!" << endl;
            const int Nl = dimensionsOfCubesToReceive(i, 0), Nc = dimensionsOfCubesToReceive(i, 1);
            recv_buff(i).resize(Nl, Nc);
            // we post the receives only if my_id == recv_proc_number, in order to not post
            // all the receives at the same time, which could deadlock
            // because of too many communications
            int BUF_SIZE = recv_buff(i).size();
            ierror = MPI_Irecv(recv_buff(i).data(), BUF_SIZE, MPI::COMPLEX, src, cubeNumber, MPI::COMM_WORLD, &Array_irecv_request(i));
          }
        }

        //if (my_id == recv_proc_number) cout << "Receiving Process " << recv_proc_number << ", Sending Process " << send_proc_number << ": All commands posted. Waiting to achieve......" << endl;
        // we wait for all the communications to finish
        if (my_id == recv_proc_number) {
          for (int i=0 ; i<NCubesToReceive ; ++i) ierror = MPI_Wait(&Array_irecv_request(i), &Array_irecv_status(i));
        }
        if (my_id == send_proc_number) {
          for (int i=0 ; i<NCubesToSend ; ++i) ierror = MPI_Wait(&Array_isend_request(i), &Array_isend_status(i));
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
  ierror = MPI_Barrier(MPI::COMM_WORLD);

  // Get peak memory usage of each rank
  long memusage_local = MemoryUsageGetPeak();
  std::cout << "MEMINFO " << argv[0] << " rank " << my_id << " mem=" << memusage_local/(1024*1024) << " MB" << std::endl;

  MPI::Finalize();
  return 0;
}

