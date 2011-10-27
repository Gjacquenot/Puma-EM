#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <random/uniform.h>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <time.h>

using namespace blitz;

const complex<double> I(0.0, 1.0);

int main(void) {
  int N = 4, ierror;
  int dest;
  int i;
  int master = 0;
  int my_id, tag = 99;
  int num_procs;
  MPI_Status status;

  MPI::Init();

  num_procs = MPI::COMM_WORLD.Get_size ( );
//
//  Get the individual process ID.
//
  my_id = MPI::COMM_WORLD.Get_rank ( );
//
//  Print a message.
//
  if ( my_id == master ) {
    cout << "Hi, I am the master node." << endl;
    cout << "###########################################" << endl;
    cout << "number of processes = " << num_procs << endl;
    cout << "###########################################" << endl;

    complex<double> a(4.0, -3.0);
    cout << "on process " << my_id << ", a = " << a << endl;
    MPI_Send(&a, 1, MPI::DOUBLE_COMPLEX, 1, tag, MPI_COMM_WORLD);

    // now a more complicated dataset
    Array<complex<double>, 1> A(N);
    ranlib::Uniform<double> a_i;
    a_i.seed((unsigned int)time(0));
    for (int i=0 ; i<A.size(); ++i) A(i) = a_i.random() + a_i.random() * I;
    cout << "on process " << my_id << ", A = " << A << endl;
    MPI_Send(A.data(), N, MPI::DOUBLE_COMPLEX, 1, tag, MPI_COMM_WORLD);

    // even more complicated
    Array<complex<float>, 2> AA(N, 2);
    for (int i=0 ; i<AA.extent(0); ++i) {
      for (int j=0 ; j<AA.extent(1); ++j) AA(i, j) = a_i.random() + a_i.random() * I;
    }
    cout << "on process " << my_id << ", AA = " << AA << endl;
    MPI_Send(AA.data(), AA.size(), MPI::COMPLEX, 1, tag, MPI_COMM_WORLD);

    // Ultimate test: use of packing
/*    Array<complex<float>, 2> AAA(N+3, 2);
    for (int i=0 ; i<AAA.extent(0); ++i) {
      for (int j=0 ; j<AAA.extent(1); ++j) AAA(i, j) = a_i.random() + a_i.random() * I;
    }
    Array<complex<float>, 1> AAABuffer(AAA.size() + 1);
    int S = AAA.size();
    int position = 0;
    MPI_Pack(&S, 1, MPI::INT, AAABuffer.data(), AAABuffer.size(), &position, MPI_COMM_WORLD);
    MPI_Pack(AAA.data(), S, MPI::COMPLEX, AAABuffer.data(), AAABuffer.size(), &position, MPI_COMM_WORLD);*/
  }
  else {
    cout << "Hi, I am node number " << my_id << endl;
    complex<double> b;
    MPI_Recv(&b, 1, MPI::DOUBLE_COMPLEX, 0, tag, MPI_COMM_WORLD, &status);
    cout << "on process " << my_id << ", b = 2.0*a = " << 2.0*b << endl;

    // now a more complicated dataset
    Array<complex<double>, 1> B(N);
    MPI_Recv(B.data(), N, MPI::DOUBLE_COMPLEX, 0, tag, MPI_COMM_WORLD, &status);
    cout << "on process " << my_id << ", B = " << B << endl;

    // even more complicated
    Array<complex<float>, 2> BB(N, 2);
    MPI_Recv(BB.data(), BB.size(), MPI::COMPLEX, 0, tag, MPI_COMM_WORLD, &status);
    cout << "on process " << my_id << ", BB = " << BB << endl;
  }

/*  N = 9;
  Array<complex<float>, 1> A(N), B, C(N);
  ranlib::Uniform<double> a_i;
  a_i.seed((unsigned int)time(0));
  for (int i=0 ; i<A.size(); ++i) A(i) = a_i.random() + a_i.random() * I;
  if ( my_id == master ) cout << "A = " << A << endl;

  Array<int, 1> MPI_Scatterv_scounts(num_procs), MPI_Scatterv_displs(num_procs);
  int displacement = 0;
  for (int i=0 ; i<num_procs ; ++i) {
    if (i<num_procs-1) MPI_Scatterv_scounts(i) = N/num_procs;
    else MPI_Scatterv_scounts(i) = N - displacement;
    MPI_Scatterv_displs(i) = displacement;
    displacement += MPI_Scatterv_scounts(i);
    if (my_id==i) B.resize(MPI_Scatterv_scounts(i));
  }
  if ( my_id == master ) cout << MPI_Scatterv_scounts << endl;
  if ( my_id == master ) cout << MPI_Scatterv_displs << endl;

  ierror = MPI_Scatterv( A.data(), MPI_Scatterv_scounts.data(), MPI_Scatterv_displs.data(), MPI::COMPLEX, B.data(), B.size(), MPI::COMPLEX, 0,  MPI_COMM_WORLD);
  cout << "on process " << my_id << ", B = " << B << endl;

  ierror = MPI_Gatherv( B.data(), B.size(), MPI::COMPLEX, C.data(), MPI_Scatterv_scounts.data(), MPI_Scatterv_displs.data(), MPI::COMPLEX, 0,  MPI_COMM_WORLD);

  cout << "on process " << my_id << ", C = " << C << endl;
*/

  MPI::Finalize();
  return 0;
}





