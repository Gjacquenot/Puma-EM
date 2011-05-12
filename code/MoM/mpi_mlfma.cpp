#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace blitz;

#include "Z_sparse_MLFMA.h"
#include "octtree.h"
#include "./iterative/iterative.h"
#include "V_E_V_H.h"
#include "mesh.h"
#include "EMConstants.h"

/****************************************************************************/
/******************************* MLFMA matvec *******************************/
/****************************************************************************/
void gatherAndRedistribute(blitz::Array<std::complex<float>, 1>& x, const int my_id, const int num_procs)
{
  int ierror, N_RWG = x.size();
  blitz::Array<std::complex<float>, 1> recvBuf(N_RWG);
  ierror = MPI_Allreduce(x.data(), recvBuf.data(), N_RWG, MPI::COMPLEX, MPI::SUM, MPI::COMM_WORLD);
  x = recvBuf;
}

class MatvecMLFMA {

  public:
    Octtree * pOcttree;
    int procNumber;
    int totalProcNumber;
    int N_RWG;
    blitz::Array<int, 1> localRWGnumbers;

    // constructors
    MatvecMLFMA(void){;};
    MatvecMLFMA(Octtree & /*octtree*/,
                const int /*numberOfRWG*/,
                const blitz::Array<int, 1>& /*localRWGindexes*/);
    // destructor
    ~MatvecMLFMA(void){localRWGnumbers.free();};
    // copy operators
    void copyMatvecMLFMA (const MatvecMLFMA&);
    MatvecMLFMA(const MatvecMLFMA&); // copy constructor
    MatvecMLFMA& operator=(const MatvecMLFMA&); // copy assignment operator

    // function
    void setProcNumber(const int n) {procNumber = n;};
    const int getProcNumber(void) const {return procNumber;};
    void setTotalProcNumber(const int n) {totalProcNumber = n;};
    const int getTotalProcNumber(void) const {return totalProcNumber;};
    const int getN_RWG(void) const {return N_RWG;};
    void matvecZnear(blitz::Array<std::complex<float>, 1> & /*y*/,
                     const blitz::Array<std::complex<float>, 1> & /*x*/);
    blitz::Array<std::complex<float>, 1> matvec(const blitz::Array<std::complex<float>, 1> & /*x*/);
};

MatvecMLFMA::MatvecMLFMA(Octtree & octtree,
                         const int numberOfRWG,
                         const blitz::Array<int, 1>& localRWGindexes)
{
  pOcttree = &octtree;
  procNumber = MPI::COMM_WORLD.Get_rank();
  totalProcNumber = MPI::COMM_WORLD.Get_size();
  N_RWG = numberOfRWG;
  localRWGnumbers.resize(localRWGindexes.size());
  localRWGnumbers = localRWGindexes;
}

void MatvecMLFMA::copyMatvecMLFMA(const MatvecMLFMA& matvecMLFMAtoCopy) // copy member function
{
  pOcttree = matvecMLFMAtoCopy.pOcttree;
  procNumber = matvecMLFMAtoCopy.procNumber;
  totalProcNumber = matvecMLFMAtoCopy.totalProcNumber;
  N_RWG = matvecMLFMAtoCopy.N_RWG;
  localRWGnumbers.resize(matvecMLFMAtoCopy.localRWGnumbers.size());
  localRWGnumbers = matvecMLFMAtoCopy.localRWGnumbers;
}

MatvecMLFMA::MatvecMLFMA(const MatvecMLFMA& matvecMLFMAtoCopy) // copy constructor
{
  copyMatvecMLFMA(matvecMLFMAtoCopy);
}

MatvecMLFMA& MatvecMLFMA::operator=(const MatvecMLFMA& matvecMLFMAtoCopy) { // copy assignment
  copyMatvecMLFMA(matvecMLFMAtoCopy);
  return *this;
}

void MatvecMLFMA::matvecZnear(blitz::Array<std::complex<float>, 1> & y, const blitz::Array<std::complex<float>, 1> & x)
{
  const int my_id = MPI::COMM_WORLD.Get_rank();
  const string pathToReadFrom = "./tmp" + intToString(my_id) + "/Z_near/", Z_name = "Z_CFIE_near";
  Array<int, 1> chunkNumbers;
  readIntBlitzArray1DFromASCIIFile(pathToReadFrom + "chunkNumbers.txt", chunkNumbers);
  Z_sparse_MLFMA Z_near;
  for (int i=0 ; i<chunkNumbers.size() ; i++) {
    int number = chunkNumbers(i);
    Z_near.setZ_sparse_MLFMAFromFile(pathToReadFrom, Z_name, number);
    Z_near.matvec_Z_PQ_near(y, x);
  }
}

blitz::Array<std::complex<float>, 1> MatvecMLFMA::matvec(const blitz::Array<std::complex<float>, 1> & x)
{
  // distribution of x among processes
  blitz::Array<std::complex<float>, 1> xTmp(this->N_RWG);
  xTmp = 0.0;
  for (int i=0 ; i<this->localRWGnumbers.size() ; ++i) xTmp(this->localRWGnumbers(i)) = x(i);
  // we should now gather and redistribute the result among the processes
  gatherAndRedistribute(xTmp, getProcNumber(), getTotalProcNumber());

  blitz::Array<std::complex<float>, 1> y(this->N_RWG);
  y = 0.0;
  // far-field multiplication
  pOcttree->ZIFarComputation(y, xTmp);
  // near-field multiplication
  matvecZnear(y, xTmp);
  xTmp.free();
  // we should now gather and redistribute the result among the processes
  gatherAndRedistribute(y, getProcNumber(), getTotalProcNumber());
  // we now select only the elements to return
  blitz::Array<std::complex<float>, 1> yTmp(this->localRWGnumbers.size());
  for (int i=0 ; i<this->localRWGnumbers.size() ; ++i) yTmp(i) = y(this->localRWGnumbers(i));
  return yTmp;
}

/****************************************************************************/
/******************************* Left Frob precond **************************/
/****************************************************************************/
class LeftFrobPsolveMLFMA {

  public:
    int procNumber;
    int totalProcNumber;
    int N_RWG;
    blitz::Array<int, 1> localRWGnumbers;

    // constructors
    LeftFrobPsolveMLFMA(void){;};
    LeftFrobPsolveMLFMA(const int /*numberOfRWG*/,
                        const blitz::Array<int, 1>& /*localRWGindexes*/);
    // destructor
    ~LeftFrobPsolveMLFMA(void){localRWGnumbers.free();};
    // copy operators
    void copyLeftFrobPsolveMLFMA (const LeftFrobPsolveMLFMA&);
    LeftFrobPsolveMLFMA(const LeftFrobPsolveMLFMA&); // copy constructor
    LeftFrobPsolveMLFMA& operator=(const LeftFrobPsolveMLFMA&); // copy assignment operator

    // function
    blitz::Array<std::complex<float>, 1> psolve(const blitz::Array<std::complex<float>, 1> & /*x*/);
};

LeftFrobPsolveMLFMA::LeftFrobPsolveMLFMA(const int numberOfRWG,
                                         const blitz::Array<int, 1>& localRWGindexes)
{
  procNumber = MPI::COMM_WORLD.Get_rank();
  totalProcNumber = MPI::COMM_WORLD.Get_size();
  N_RWG = numberOfRWG;
  localRWGnumbers.resize(localRWGindexes.size());
  localRWGnumbers = localRWGindexes;
}

void LeftFrobPsolveMLFMA::copyLeftFrobPsolveMLFMA(const LeftFrobPsolveMLFMA& leftFrobPsolveMLFMAtoCopy) // copy member function
{
  procNumber = leftFrobPsolveMLFMAtoCopy.procNumber;
  totalProcNumber = leftFrobPsolveMLFMAtoCopy.totalProcNumber;
  N_RWG = leftFrobPsolveMLFMAtoCopy.N_RWG;
  localRWGnumbers.resize(leftFrobPsolveMLFMAtoCopy.localRWGnumbers.size());
  localRWGnumbers = leftFrobPsolveMLFMAtoCopy.localRWGnumbers;
}

LeftFrobPsolveMLFMA::LeftFrobPsolveMLFMA(const LeftFrobPsolveMLFMA& leftFrobPsolveMLFMAtoCopy) // copy constructor
{
  copyLeftFrobPsolveMLFMA(leftFrobPsolveMLFMAtoCopy);
}

LeftFrobPsolveMLFMA& LeftFrobPsolveMLFMA::operator=(const LeftFrobPsolveMLFMA& leftFrobPsolveMLFMAtoCopy) { // copy assignment
  copyLeftFrobPsolveMLFMA(leftFrobPsolveMLFMAtoCopy);
  return *this;
}

blitz::Array<std::complex<float>, 1> LeftFrobPsolveMLFMA::psolve(const blitz::Array<std::complex<float>, 1> & x)
{
  const int my_id = MPI::COMM_WORLD.Get_rank();
  blitz::Array<std::complex<float>, 1> xTmp(this->N_RWG);
  xTmp = 0.0;
  for (int i=0 ; i<localRWGnumbers.size() ; ++i) xTmp(localRWGnumbers(i)) = x(i);
  gatherAndRedistribute(xTmp, procNumber, totalProcNumber);

  blitz::Array<std::complex<float>, 1> y(this->N_RWG);
  y = 0.0;
  const string pathToReadFrom = "./tmp" + intToString(my_id) + "/Mg_LeftFrob/", Z_name = "Mg_LeftFrob";
  blitz::Array<int, 1> chunkNumbers;
  readIntBlitzArray1DFromASCIIFile(pathToReadFrom + "chunkNumbers.txt", chunkNumbers);
  Z_sparse_MLFMA Mg_LeftFrob;
  for (int i=0 ; i<chunkNumbers.size() ; i++) {
    int number = chunkNumbers(i);
    Mg_LeftFrob.setZ_sparse_MLFMAFromFile(pathToReadFrom, Z_name, number);
    //Mg_LeftFrob.printZ_CFIE_near();
    Mg_LeftFrob.matvec_Z_PQ_near(y, xTmp);
  }
  xTmp.free();
  // we should now gather and redistribute the result among the processes
  gatherAndRedistribute(y, procNumber, totalProcNumber);
  // we now select only the elements to return
  blitz::Array<std::complex<float>, 1> yTmp(localRWGnumbers.size());
  for (int i=0 ; i<localRWGnumbers.size() ; ++i) yTmp(i) = y(localRWGnumbers(i));
  return yTmp;
}

/****************************************************************************/
/*************************** AMLFMA Preconditioner **************************/
/****************************************************************************/
class PsolveAMLFMA { // should use FGMRES and _NOT_ GMRES or BiCGSTAB
    MatvecMLFMA *pMatvecMLFMA;
    LeftFrobPsolveMLFMA *pLeftFrobPsolveMLFMA;
    int inner_maxiter;
    int inner_restart;
    int my_id;
    int num_procs;
    int N_RWG;
    int NlocalRWG;
    float inner_tol;
    string inner_solver;

  public:
    // constructors
    PsolveAMLFMA(MatvecMLFMA & /*matvecMLFMA*/,
                 LeftFrobPsolveMLFMA & /*leftFrobPsolveMLFMA*/,
                 const float /*INNER_TOL*/,
                 const int /*INNER_MAXITER*/,
                 const int /*INNER_RESTART*/,
                 const int /*numberOfRWG*/,
                 const string & /*INNER_SOLVER*/);
    ~PsolveAMLFMA(void){};
    // function
    blitz::Array<std::complex<float>, 1> psolve(const blitz::Array<std::complex<float>, 1> & /*x*/);
};

PsolveAMLFMA::PsolveAMLFMA (MatvecMLFMA & matvecMLFMA,
                            LeftFrobPsolveMLFMA & leftFrobPsolveMLFMA,
                            const float INNER_TOL,
                            const int INNER_MAXITER,
                            const int INNER_RESTART,
                            const int numberOfRWG,
                            const string & INNER_SOLVER)
{
  pMatvecMLFMA = &matvecMLFMA;
  pLeftFrobPsolveMLFMA = &leftFrobPsolveMLFMA;
  inner_maxiter = INNER_MAXITER;
  inner_restart = INNER_RESTART;
  inner_tol = INNER_TOL;
  my_id = MPI::COMM_WORLD.Get_rank();
  num_procs = MPI::COMM_WORLD.Get_size();
  N_RWG = numberOfRWG;
  NlocalRWG = matvecMLFMA.localRWGnumbers.size();
  inner_solver = INNER_SOLVER;
}

blitz::Array<std::complex<float>, 1> PsolveAMLFMA::psolve(const blitz::Array<std::complex<float>, 1>& x)
{
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(pMatvecMLFMA, &MatvecMLFMA::matvec);
  PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > innerPsolve(pLeftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);

  const string TMP = "./tmp" + intToString(my_id), ITERATIVE_DATA_PATH = TMP + "/iterative_data/";
  int iter, flag;
  double error;
  blitz::Array<std::complex<float>, 1> y(this->NlocalRWG);
  y = 0.0;
  if (inner_solver == "BICGSTAB") bicgstab(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, 100, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  else if (inner_solver == "GMRES") gmres(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, inner_restart, inner_maxiter, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  else if (inner_solver == "RGMRES") fgmres(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, inner_restart, inner_maxiter, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  else if (inner_solver == "FGMRES") fgmres(y, error, iter, flag, matvec, innerPsolve, x, inner_tol, inner_restart, inner_maxiter, my_id, num_procs, ITERATIVE_DATA_PATH + "inner_convergence.txt");
  return y;
}

/****************************************************************************/
/************************* some useful functions ****************************/
/****************************************************************************/

void computeE_obs(blitz::Array<std::complex<double>, 2>& E_obs,
                  const blitz::Array<double, 2>& r_obs,
                  const LocalMesh &  local_target_mesh,
                  const blitz::Array<std::complex<float>, 1>& ZI,
                  const std::complex<float> eps_r,
                  const std::complex<float> mu_r,
                  const float w)
{
  int ierror, FULL_PRECISION = 1;
  blitz::Range all = blitz::Range::all();
  blitz::Array<std::complex<float>, 1> CFIE_for_tE(4);
  CFIE_for_tE = 1.0, 0.0, 0.0, 0.0;
  E_obs.resize(r_obs.extent(0), r_obs.extent(1));
  for (int j=0 ; j<r_obs.extent(0) ; ++j) {
    for (int i=0 ; i<3 ; ++i) {
      blitz::Array<std::complex<float>, 1> V_tE;
      blitz::Array<std::complex<double>, 1> J_obs(3);
      J_obs = 0.0;
      J_obs(i) = 1.0;
      local_V_CFIE_dipole (V_tE, J_obs, r_obs(j, all), local_target_mesh, w, eps_r, mu_r, CFIE_for_tE, FULL_PRECISION);
      std::complex<float> local_e_tmp = sum(V_tE * ZI), e_tmp;
      ierror = MPI_Allreduce(&local_e_tmp, &e_tmp, 1, MPI::COMPLEX, MPI::SUM, MPI::COMM_WORLD);
      E_obs(j, i) = e_tmp;
    }
  }
}

void computeForOneExcitation(Octtree & octtree,
                             LocalMesh & local_target_mesh,
                             const string SOLVER,
                             const string TMP,
                             const string OCTTREE_DATA_PATH,
                             const string MESH_DATA_PATH,
                             const string V_CFIE_DATA_PATH,
                             const string RESULT_DATA_PATH,
                             const string ITERATIVE_DATA_PATH)
{
  blitz::Range all = blitz::Range::all();
  int num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0, N_local_RWG = local_target_mesh.N_local_RWG;
  const float w = octtree.w;
  const std::complex<float> eps_r = octtree.eps_r, mu_r = octtree.mu_r;
  blitz::Array<int, 1> localRWGNumbers(local_target_mesh.localRWGNumbers.size());
  localRWGNumbers = local_target_mesh.localRWGNumbers;
  // functors declarations
  const string PRECOND = "FROB";
  int N_RWG, iter, flag, RESTART, MAXITER;
  double error, TOL;
  readIntFromASCIIFile(OCTTREE_DATA_PATH + "N_RWG.txt", N_RWG);
  MatvecMLFMA matvecMLFMA(octtree, N_RWG, localRWGNumbers);
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(&matvecMLFMA, &MatvecMLFMA::matvec);
  LeftFrobPsolveMLFMA leftFrobPsolveMLFMA(N_RWG, localRWGNumbers);

  // excitation : V_CFIE computation
  blitz::Array<std::complex<float>, 1> V_CFIE(N_local_RWG);
  V_CFIE = 0.0;
  int DIPOLES_EXCITATION, PLANE_WAVE_EXCITATION, V_FULL_PRECISION;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "DIPOLES_EXCITATION.txt", DIPOLES_EXCITATION);
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "PLANE_WAVE_EXCITATION.txt", PLANE_WAVE_EXCITATION);
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "V_FULL_PRECISION.txt", V_FULL_PRECISION);
  if (DIPOLES_EXCITATION==1) {
    blitz::Array<std::complex<double>, 2> J_dip, M_dip;
    blitz::Array<double, 2> r_J_dip, r_M_dip;
    int J_DIPOLES_EXCITATION, M_DIPOLES_EXCITATION;
    // electric dipoles
    readIntFromASCIIFile(V_CFIE_DATA_PATH + "J_DIPOLES_EXCITATION.txt", J_DIPOLES_EXCITATION);
    if (J_DIPOLES_EXCITATION==1) {
      readBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "J_dip.txt", J_dip);
      readBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "r_J_dip.txt", r_J_dip);
      if (my_id==0) {
        //cout << "J_dip.txt = " << J_dip << endl;
        //cout << "r_J_dip.txt = " << r_J_dip << endl;
      }
      blitz::Array<std::complex<float>, 1> V_CFIE_tmp;
      const char CURRENT_TYPE = 'J';
      local_V_CFIE_dipole_array (V_CFIE_tmp, J_dip, r_J_dip, local_target_mesh, w, eps_r, mu_r, octtree.CFIE, CURRENT_TYPE, V_FULL_PRECISION);
      V_CFIE += V_CFIE_tmp;
    }
    // magnetic dipoles
    readIntFromASCIIFile(V_CFIE_DATA_PATH + "M_DIPOLES_EXCITATION.txt", M_DIPOLES_EXCITATION);
    if (M_DIPOLES_EXCITATION==1) {
      readBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "M_dip.txt", M_dip);
      readBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "r_M_dip.txt", r_M_dip);
      if (my_id==0) {
        //cout << "M_dip.txt = " << M_dip << endl;
        //cout << "r_M_dip.txt = " << r_M_dip << endl;
      }
      blitz::Array<std::complex<float>, 1> V_CFIE_tmp;
      const char CURRENT_TYPE = 'M';
      local_V_CFIE_dipole_array (V_CFIE_tmp, J_dip, r_J_dip, local_target_mesh, w, eps_r, mu_r, octtree.CFIE, CURRENT_TYPE, V_FULL_PRECISION);
      V_CFIE += V_CFIE_tmp;
    }
  }
  if (PLANE_WAVE_EXCITATION==1) {
    double theta_inc, phi_inc;
    readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "theta_inc.txt", theta_inc);
    readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "phi_inc.txt", phi_inc);
    // local coordinate system
    blitz::Array<double, 1> r_hat(3), theta_hat(3), phi_hat(3), r_ref(3);
    r_hat = sin(theta_inc)*cos(phi_inc), sin(theta_inc)*sin(phi_inc), cos(theta_inc);
    theta_hat = cos(theta_inc)*cos(phi_inc), cos(theta_inc)*sin(phi_inc), -sin(theta_inc);
    phi_hat = -sin(phi_inc), cos(phi_inc), 0.0;
    blitz::Array<double, 1> k_hat(-1.0 * r_hat);
    r_ref = 0.0;
    // excitation electric field
    // we now read the incoming field amplitude and phase and polarization
    // only 2 components are needed: E_theta and E_phi
    blitz::Array<std::complex<double>, 1> E_inc_components(2), E_inc(3);
    readComplexDoubleBlitzArray1DFromASCIIFile( V_CFIE_DATA_PATH + "E_inc.txt", E_inc_components );
    E_inc = E_inc_components(0) * theta_hat + E_inc_components(1) * phi_hat;
    blitz::Array<std::complex<float>, 1> V_CFIE_tmp;
    local_V_CFIE_plane (V_CFIE_tmp, E_inc, k_hat, r_ref, local_target_mesh, octtree.w, octtree.eps_r, octtree.mu_r, octtree.CFIE, V_FULL_PRECISION);
    V_CFIE += V_CFIE_tmp;
  }
  local_target_mesh.resizeToZero();
  // iterative solving
  readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "TOL.txt", TOL);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "RESTART.txt", RESTART);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "MAXITER.txt", MAXITER);
  blitz::Array<std::complex<float>, 1> ZI(N_local_RWG);
  ZI = 0.0;
  if (SOLVER=="BICGSTAB") {
    PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
    bicgstab(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  }
  else if (SOLVER=="GMRES") {
    PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
    gmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  } 
  else if (SOLVER=="RGMRES") {
    PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
    fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  } 
  else if (SOLVER=="FGMRES") {
    string INNER_SOLVER;
    readStringFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_SOLVER.txt", INNER_SOLVER);
    int INNER_MAXITER, INNER_RESTART;
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_MAXITER.txt", INNER_MAXITER);
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_RESTART.txt", INNER_RESTART);
    double INNER_TOL;
    readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_TOL.txt", INNER_TOL);
    PsolveAMLFMA psolveAMLFMA(matvecMLFMA, leftFrobPsolveMLFMA, INNER_TOL, INNER_MAXITER, INNER_RESTART, N_RWG, INNER_SOLVER);
    PrecondFunctor< std::complex<float>, PsolveAMLFMA > psolve(&psolveAMLFMA, &PsolveAMLFMA::psolve);
    fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
  }
  else {
    cout << "Bad solver choice!! Solver is BICGSTAB or (F)GMRES, and you chose " << SOLVER << endl;
    exit(1);
  }
  octtree.resizeSdownLevelsToZero();

  // now computing the E_field at the r_obs
  blitz::Array<double, 2> r_obs;
  readDoubleBlitzArray2DFromASCIIFile( V_CFIE_DATA_PATH + "r_obs.txt", r_obs);
  blitz::Array<std::complex<double>, 2> E_obs;
  local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
  computeE_obs(E_obs, r_obs, local_target_mesh, ZI, eps_r, mu_r, w);
  local_target_mesh.resizeToZero();
  if (my_id==master) writeComplexDoubleBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "E_obs.txt", E_obs);
  if (my_id==master) writeDoubleBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "r_obs.txt", r_obs);

  // we now gather all the ZI from the other processes
  blitz::Array<std::complex<float>, 1> ZI_tmp(N_RWG);
  ZI_tmp = 0.0;
  for (int i=0 ; i<N_local_RWG ; ++i) ZI_tmp(localRWGNumbers(i)) = ZI(i);
  gatherAndRedistribute(ZI_tmp, my_id, num_procs);
  // calculating the far fields
  blitz::Array<float, 1> octtreeXthetas_coarsest, octtreeXphis_coarsest;
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXphis_coarsest.txt", octtreeXphis_coarsest);
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXthetas_coarsest.txt", octtreeXthetas_coarsest);
  blitz::Array<std::complex<float>, 2> e_theta_far, e_phi_far;
  octtree.computeFarField(e_theta_far, e_phi_far, octtreeXthetas_coarsest, octtreeXphis_coarsest, ZI_tmp, OCTTREE_DATA_PATH);
  if (my_id==master) {
    writeComplexFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "e_theta_far_ASCII.txt", e_theta_far);
    writeComplexFloatBlitzArray2DToBinaryFile(RESULT_DATA_PATH + "e_theta_far_Binary.txt", e_theta_far);
    writeComplexFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "e_phi_far_ASCII.txt", e_phi_far);
    writeComplexFloatBlitzArray2DToBinaryFile(RESULT_DATA_PATH + "e_phi_far_Binary.txt", e_phi_far);
    // thetas, phis
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "phis_far_field_ASCII.txt", octtreeXphis_coarsest);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "thetas_far_field_ASCII.txt", octtreeXthetas_coarsest);
  }
  // we now write the solution to disk...
  if ( my_id == master ) {
    string filename = TMP + "/ZI/ZI.txt";
    ofstream ofs(filename.c_str(), blitz::ios::binary);
    ofs.write((char *)(ZI_tmp.data()), ZI_tmp.size()*8);
    ofs.close();
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "numberOfMatvecs.txt", octtree.getNumberOfUpdates());
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "iter.txt", iter);
    cout << endl;
  }
}

void computeMonostaticRCS(Octtree & octtree,
                          LocalMesh & local_target_mesh,
                          const string SOLVER,
                          const string TMP,
                          const string OCTTREE_DATA_PATH,
                          const string MESH_DATA_PATH,
                          const string V_CFIE_DATA_PATH,
                          const string RESULT_DATA_PATH,
                          const string ITERATIVE_DATA_PATH)
{
  blitz::Range all = blitz::Range::all();
  int num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0, N_local_RWG = local_target_mesh.N_local_RWG;
  
  const std::complex<float> eps_r = octtree.eps_r, mu_r = octtree.mu_r;
  blitz::Array<int, 1> localRWGNumbers(local_target_mesh.localRWGNumbers.size());
  localRWGNumbers = local_target_mesh.localRWGNumbers;

  // functors declarations
  const string PRECOND = "FROB";
  int N_RWG, iter, flag, RESTART, MAXITER, USE_PREVIOUS_SOLUTION, MONOSTATIC_BY_BISTATIC_APPROX, V_FULL_PRECISION;
  double error, TOL, MAX_DELTA_PHASE;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "V_FULL_PRECISION.txt", V_FULL_PRECISION);
  readIntFromASCIIFile(OCTTREE_DATA_PATH + "N_RWG.txt", N_RWG);
  readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "TOL.txt", TOL);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "RESTART.txt", RESTART);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "MAXITER.txt", MAXITER);
  readIntFromASCIIFile(TMP + "/USE_PREVIOUS_SOLUTION.txt", USE_PREVIOUS_SOLUTION);
  readIntFromASCIIFile(TMP + "/MONOSTATIC_BY_BISTATIC_APPROX.txt", MONOSTATIC_BY_BISTATIC_APPROX);
  readDoubleFromASCIIFile(TMP + "/MAXIMUM_DELTA_PHASE.txt", MAX_DELTA_PHASE);
  MAX_DELTA_PHASE *= M_PI/180.0;
  if (MONOSTATIC_BY_BISTATIC_APPROX!=1) MAX_DELTA_PHASE = 0.0;

  MatvecMLFMA matvecMLFMA(octtree, N_RWG, localRWGNumbers);
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(&matvecMLFMA, &MatvecMLFMA::matvec);
  LeftFrobPsolveMLFMA leftFrobPsolveMLFMA(N_RWG, localRWGNumbers);
  PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  if (SOLVER=="FGMRES") {
    string INNER_SOLVER;
    readStringFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_SOLVER.txt", INNER_SOLVER);
    int INNER_MAXITER, INNER_RESTART;
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_MAXITER.txt", INNER_MAXITER);
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_RESTART.txt", INNER_RESTART);
    double INNER_TOL;
    readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_TOL.txt", INNER_TOL);
    PsolveAMLFMA psolveAMLFMA(matvecMLFMA, leftFrobPsolveMLFMA, INNER_TOL, INNER_MAXITER, INNER_RESTART, N_RWG, INNER_SOLVER);
    PrecondFunctor< std::complex<float>, PsolveAMLFMA > psolve(&psolveAMLFMA, &PsolveAMLFMA::psolve);
  }
  else PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  // getting the angles at which monostatic RCS must be computed
  blitz::Array<float, 1> octtreeXthetas_coarsest, octtreeXphis_coarsest;
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXphis_coarsest.txt", octtreeXphis_coarsest);
  readFloatBlitzArray1DFromASCIIFile(OCTTREE_DATA_PATH + "octtreeXthetas_coarsest.txt", octtreeXthetas_coarsest);
  const int N_theta(octtreeXthetas_coarsest.size()), N_phi(octtreeXphis_coarsest.size());
  blitz::Array<float, 2> RCS_VV(N_theta, N_phi), RCS_HH(N_theta, N_phi), RCS_HV(N_theta, N_phi);
  RCS_VV = 1.0;
  RCS_HH = 1.0;
  RCS_HV = 1.0;
  int COMPUTE_RCS_HH, COMPUTE_RCS_HV, COMPUTE_RCS_VV;
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HH.txt", COMPUTE_RCS_HH);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HV.txt", COMPUTE_RCS_HV);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_VV.txt", COMPUTE_RCS_VV);
  // r_ref
  blitz::Array<double, 1> r_ref(3);
  for (int i=0 ; i<3 ; ++i) r_ref(i) = octtree.big_cube_center_coord(i);
  // we now make the calculations of the parameters for the MONOSTATIC_BY_BISTATIC_APPROX
  // see pages 94-96 of Chew 2001, "Fast and Efficient Methods in Computational EM"
  const double big_cube_side_length = 2.0 * (octtree.big_cube_center_coord(0) - octtree.big_cube_lower_coord(0));
  const double LL = big_cube_side_length * sqrt(2.0);
  const double Delta_Phi = octtreeXphis_coarsest(1) - octtreeXphis_coarsest(0);
  double Beta = sqrt(4.0 * c * MAX_DELTA_PHASE / (LL * octtree.w));
  // number of phi points encompassed by Beta: number of points for which I can approximate
  // monostatic RCS by bistatic RCS
  const int BetaPoints = static_cast<int>(floor(Beta/Delta_Phi)) + 1;
  Beta = (BetaPoints-1) * Delta_Phi;
  if (my_id==master) cout << "number of BetaPoints for monostatic-bistatic approximation = " << BetaPoints << endl;
  // loop for monostatic sigma computation
  for (int excitation=0 ; excitation<2 ; ++excitation) { // 0 for H, 1 for V
    const bool HH = ((excitation==0) && (COMPUTE_RCS_HH==1));
    const bool HV = ((excitation==0) && (COMPUTE_RCS_HV==1));
    const bool VV = ((excitation==1) && (COMPUTE_RCS_VV==1));
    const bool cond = (HH || HV || VV);
    if (cond) {
      for (int t=0 ; t<N_theta ; ++t) {
        blitz::Array<std::complex<float>, 1> ZI(N_local_RWG);
        ZI = 0.0;
        const float theta = octtreeXthetas_coarsest(t);
        float phi_inc = octtreeXphis_coarsest(0) + Beta/2.0;
        int startIndexPhi = 0;
        while (startIndexPhi<N_phi) {
          octtree.resizeSdownLevelsToZero();
          if (my_id==master) {
            if (HH || HV) cout << "\nHH and HV, theta = "<< theta * 180.0/M_PI << ", phi = " << phi_inc * 180.0/M_PI << endl;
            else cout << "\nVV, theta = "<< theta * 180.0/M_PI << ", phi = " << phi_inc * 180.0/M_PI << endl;
            flush(cout);
          }
          // local coordinate system
          blitz::Array<double, 1> r_hat(3), theta_hat(3), phi_hat(3);
          r_hat = sin(theta)*cos(phi_inc), sin(theta)*sin(phi_inc), cos(theta);
          theta_hat = cos(theta)*cos(phi_inc), cos(theta)*sin(phi_inc), -sin(theta);
          phi_hat = -sin(phi_inc), cos(phi_inc), 0.0;
          blitz::Array<double, 1> k_hat(-1.0 * r_hat);

          // excitation field
          blitz::Array<std::complex<double>, 1> E_0(3);
          if (HH || HV) E_0 = 100.0 * (phi_hat + I * 0.0);
          else E_0 = -100.0 * (theta_hat + I * 0.0);
          blitz::Array<std::complex<float>, 1> V_CFIE;
          local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
          local_V_CFIE_plane (V_CFIE, E_0, k_hat, r_ref, local_target_mesh, octtree.w, octtree.eps_r, octtree.mu_r, octtree.CFIE, V_FULL_PRECISION);
          local_target_mesh.resizeToZero();
          // solving
          octtree.setNumberOfUpdates(0);
          if (USE_PREVIOUS_SOLUTION != 1) ZI = 0.0;
          if (SOLVER=="BICGSTAB") bicgstab(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if (SOLVER=="GMRES") gmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if ((SOLVER=="RGMRES") || (SOLVER=="FGMRES")) fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else {
            cout << "Bad solver choice!! Solver is BICGSTAB or (F)GMRES, and you chose " << SOLVER << endl;
            exit(1);
          }
          // far field computation
          blitz::Array<std::complex<float>, 2> e_theta_far, e_phi_far;
          blitz::Array<float, 1> thetas(1), phis(BetaPoints);
          thetas(0) = theta;
          if ((BetaPoints%2)!=0) { // if BetaPoints is odd and greater than 2
            const float space = 2.0*Delta_Phi;
            for (int j=0 ; j<BetaPoints ; ++j) phis(j) = phi_inc - Beta + j * space;
          }
          else {
            const float space = Delta_Phi;
            for (int j=0 ; j<BetaPoints/2 ; ++j) {
              phis(j) = phi_inc - Beta + j * space;
              phis(BetaPoints-1-j) = phi_inc + Beta - j * space;
            }
          }
          blitz::Array<std::complex<float>, 1> ZI_tmp(N_RWG);
          ZI_tmp = 0.0;
          for (int i=0 ; i<N_local_RWG ; ++i) ZI_tmp(localRWGNumbers(i)) = ZI(i);
          gatherAndRedistribute(ZI_tmp, my_id, num_procs);
          octtree.computeFarField(e_theta_far, e_phi_far, thetas, phis, ZI_tmp, OCTTREE_DATA_PATH);
          // filling of the RCS Arrays
          if (HH || HV) {
            for (int j=0 ; j<BetaPoints ; ++j) {
              RCS_HH(t, startIndexPhi + j) = real(e_phi_far(0, j) * conj(e_phi_far(0, j)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
              RCS_HV(t, startIndexPhi + j) = real(e_theta_far(0, j) * conj(e_theta_far(0, j)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
            }
          }
          else {
            for (int j=0 ; j<BetaPoints ; ++j) {
              RCS_VV(t, startIndexPhi + j) = real(e_theta_far(0, j) * conj(e_theta_far(0, j)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
            }
          }
          // phi update
          startIndexPhi += BetaPoints;
          phi_inc += BetaPoints * Delta_Phi; // += Beta;
          if (my_id==master) {
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HH_ASCII.txt", RCS_HH);
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HV_ASCII.txt", RCS_HV);
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_VV_ASCII.txt", RCS_VV);
          }
        }
      }
    }
  }
  if (my_id==master) {
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HH_ASCII.txt", RCS_HH);
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_HV_ASCII.txt", RCS_HV);
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "RCS_VV_ASCII.txt", RCS_VV);
    // thetas, phis
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "phis_far_field_ASCII.txt", octtreeXphis_coarsest);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "thetas_far_field_ASCII.txt", octtreeXthetas_coarsest);
    // final writes
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "numberOfMatvecs.txt", octtree.getNumberOfUpdates());
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "iter.txt", iter);
  }
}

void computeMonostaticSAR(Octtree & octtree,
                          LocalMesh & local_target_mesh,
                          const string SOLVER,
                          const string TMP,
                          const string OCTTREE_DATA_PATH,
                          const string MESH_DATA_PATH,
                          const string V_CFIE_DATA_PATH,
                          const string RESULT_DATA_PATH,
                          const string ITERATIVE_DATA_PATH)
{
  blitz::Range all = blitz::Range::all();
  int num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();
  const int master = 0, N_local_RWG = local_target_mesh.N_local_RWG;
  const double w = octtree.w;
  blitz::Array<int, 1> localRWGNumbers(local_target_mesh.localRWGNumbers.size());
  localRWGNumbers = local_target_mesh.localRWGNumbers;
  const std::complex<double> eps_r = octtree.eps_r, mu_r = octtree.mu_r;
  const std::complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps*mu);

  // functors declarations
  const string PRECOND = "FROB";
  int N_RWG, iter, flag, RESTART, MAXITER, USE_PREVIOUS_SOLUTION, MONOSTATIC_BY_BISTATIC_APPROX, V_FULL_PRECISION;
  double error, TOL, MAX_DELTA_PHASE;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "V_FULL_PRECISION.txt", V_FULL_PRECISION);
  readIntFromASCIIFile(OCTTREE_DATA_PATH + "N_RWG.txt", N_RWG);
  readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "TOL.txt", TOL);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "RESTART.txt", RESTART);
  readIntFromASCIIFile(ITERATIVE_DATA_PATH + "MAXITER.txt", MAXITER);
  readIntFromASCIIFile(TMP + "/USE_PREVIOUS_SOLUTION.txt", USE_PREVIOUS_SOLUTION);

  MatvecMLFMA matvecMLFMA(octtree, N_RWG, localRWGNumbers);
  MatvecFunctor< std::complex<float>, MatvecMLFMA > matvec(&matvecMLFMA, &MatvecMLFMA::matvec);
  LeftFrobPsolveMLFMA leftFrobPsolveMLFMA(N_RWG, localRWGNumbers);
  PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  if (SOLVER=="FGMRES") {
    string INNER_SOLVER;
    readStringFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_SOLVER.txt", INNER_SOLVER);
    int INNER_MAXITER, INNER_RESTART;
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_MAXITER.txt", INNER_MAXITER);
    readIntFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_RESTART.txt", INNER_RESTART);
    double INNER_TOL;
    readDoubleFromASCIIFile(ITERATIVE_DATA_PATH + "INNER_TOL.txt", INNER_TOL);
    PsolveAMLFMA psolveAMLFMA(matvecMLFMA, leftFrobPsolveMLFMA, INNER_TOL, INNER_MAXITER, INNER_RESTART, N_RWG, INNER_SOLVER);
    PrecondFunctor< std::complex<float>, PsolveAMLFMA > psolve(&psolveAMLFMA, &PsolveAMLFMA::psolve);
  }
  else PrecondFunctor< std::complex<float>, LeftFrobPsolveMLFMA > psolve(&leftFrobPsolveMLFMA, &LeftFrobPsolveMLFMA::psolve);
  // getting the positions at which monostatic SAR must be computed
  blitz::Array<double, 1> SAR_local_x_hat(3), SAR_local_y_hat(3), SAR_plane_origin(3);
  readDoubleBlitzArray1DFromASCIIFile(V_CFIE_DATA_PATH + "SAR_local_x_hat.txt", SAR_local_x_hat);
  readDoubleBlitzArray1DFromASCIIFile(V_CFIE_DATA_PATH + "SAR_local_y_hat.txt", SAR_local_y_hat);
  // normalization
  SAR_local_x_hat = SAR_local_x_hat / sqrt(sum(SAR_local_x_hat * SAR_local_x_hat));
  SAR_local_y_hat = SAR_local_y_hat / sqrt(sum(SAR_local_y_hat * SAR_local_y_hat));
  // other SAR data
  readDoubleBlitzArray1DFromASCIIFile(V_CFIE_DATA_PATH + "SAR_plane_origin.txt", SAR_plane_origin);
  double SAR_x_span, SAR_x_span_offset, SAR_y_span, SAR_y_span_offset;
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_x_span.txt", SAR_x_span);
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_x_span_offset.txt", SAR_x_span_offset);
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_y_span.txt", SAR_y_span);
  readDoubleFromASCIIFile(V_CFIE_DATA_PATH + "SAR_y_span_offset.txt", SAR_y_span_offset);
  int SAR_N_x_points, SAR_N_y_points;
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "SAR_N_x_points.txt", SAR_N_x_points);
  readIntFromASCIIFile(V_CFIE_DATA_PATH + "SAR_N_y_points.txt", SAR_N_y_points);
  // protection against stupidity
  if (SAR_N_x_points<1) SAR_N_x_points = 1;
  if (SAR_N_y_points<1) SAR_N_y_points = 1;
  // step increments
  double Delta_x, Delta_y;
  Delta_x = (SAR_N_x_points>1) ? SAR_x_span / (SAR_N_x_points-1) : SAR_x_span/2.0;
  Delta_y = (SAR_N_y_points>1) ? SAR_y_span / (SAR_N_y_points-1) : SAR_y_span/2.0;

  blitz::Array<float, 2> r_SAR(SAR_N_x_points * SAR_N_y_points, 3);
  blitz::Array<float, 1> RCS_VV(SAR_N_x_points * SAR_N_y_points), RCS_HH(SAR_N_x_points * SAR_N_y_points), RCS_HV(SAR_N_x_points * SAR_N_y_points);
  r_SAR = 0.0;
  RCS_VV = 1.0;
  RCS_HH = 1.0;
  RCS_HV = 1.0;
  int COMPUTE_RCS_HH, COMPUTE_RCS_HV, COMPUTE_RCS_VV;
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HH.txt", COMPUTE_RCS_HH);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_HV.txt", COMPUTE_RCS_HV);
  readIntFromASCIIFile(TMP + "/COMPUTE_RCS_VV.txt", COMPUTE_RCS_VV);
  // r_ref
  blitz::Array<double, 1> r_ref(3);
  for (int i=0 ; i<3 ; ++i) r_ref(i) = octtree.big_cube_center_coord(i);
  // loop for monostatic sigma computation
  for (int excitation=0 ; excitation<2 ; ++excitation) { // 0 for H, 1 for V
    const bool HH = ((excitation==0) && (COMPUTE_RCS_HH==1));
    const bool HV = ((excitation==0) && (COMPUTE_RCS_HV==1));
    const bool VV = ((excitation==1) && (COMPUTE_RCS_VV==1));
    const bool cond = (HH || HV || VV);
    if (cond) {
      for (int t=0 ; t<SAR_N_y_points ; ++t) {
        blitz::Array<std::complex<float>, 1> ZI(N_local_RWG);
        ZI = 0.0;
        for (int s=0 ; s<SAR_N_x_points ; ++s) {
          int index = t * SAR_N_x_points + s;
          blitz::Array<double, 1> r_src(3);
          r_src = SAR_plane_origin + SAR_x_span_offset * SAR_local_x_hat + SAR_y_span_offset * SAR_local_y_hat;
          r_src += s * Delta_x * SAR_local_x_hat + t * Delta_y * SAR_local_y_hat;
          r_SAR(index, all) = 1.0 * r_src;
          octtree.resizeSdownLevelsToZero();
          if (my_id==master) {
            if (HH || HV) blitz::cout << "\nHH and HV, r_ant = "<< r_src << blitz::endl;
            else blitz::cout << "\nVV, r_ant = "<< r_src << blitz::endl;
            blitz::flush(blitz::cout);
          }
          // excitation field
          blitz::Array<std::complex<double>, 1> J_ant(3), E_0(3);
          if (HH || HV) J_ant = 100. * SAR_local_x_hat;
          else J_ant = 100. * SAR_local_y_hat;
          blitz::Array<std::complex<float>, 1> V_CFIE;
          local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
          local_V_CFIE_dipole (V_CFIE, J_ant, r_src, local_target_mesh, w, eps_r, mu_r, octtree.CFIE, V_FULL_PRECISION);
          local_target_mesh.resizeToZero();
          // target incoming field: reference field for RCS computation
          blitz::Array<std::complex<double>, 2> G_EJ(3, 3), G_HJ(3, 3);
          blitz::TinyVector<double, 3> r_dip, r_obs2;
          for (int m=0 ; m<3 ; m++) {
            r_dip(m) = r_src(m);
            r_obs2(m) = r_ref(m);
          }
          G_EJ_G_HJ (G_EJ, G_HJ, r_dip, r_obs2, eps, mu, k);
          for (int m=0 ; m<3 ; m++) E_0(m) = G_EJ(m, 0) * J_ant(0) + G_EJ(m, 1) * J_ant(1) + G_EJ(m, 2) * J_ant(2);
          // solving
          octtree.setNumberOfUpdates(0);
          if (USE_PREVIOUS_SOLUTION != 1) ZI = 0.0;
          if (SOLVER=="BICGSTAB") bicgstab(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if (SOLVER=="GMRES") gmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else if ((SOLVER=="RGMRES") || (SOLVER=="FGMRES")) fgmres(ZI, error, iter, flag, matvec, psolve, V_CFIE, TOL, RESTART, MAXITER, my_id, num_procs, ITERATIVE_DATA_PATH + "/convergence.txt");
          else {
            cout << "Bad solver choice!! Solver is BICGSTAB or (F)GMRES, and you chose " << SOLVER << endl;
            exit(1);
          }
          // field computation (O(N) version)
          blitz::Array<double, 2> r_obs(1, 3);
          blitz::Array<std::complex<double>, 2> E_obs(1, 3);
          r_obs(0, all) = r_ref;
          local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
          computeE_obs(E_obs, r_obs, local_target_mesh, ZI, static_cast<complex<float> >(eps_r), static_cast<complex<float> >(mu_r), w);
          local_target_mesh.resizeToZero();
          // filling of the RCS Arrays
          if (HH || HV) {
            blitz::Array<std::complex<double>, 1> E_H(3), E_V(3);
            E_H = sum(SAR_local_x_hat * E_obs(0, all));
            E_V = sum(SAR_local_y_hat * E_obs(0, all));
            RCS_HH(index) = real(sum(E_H * conj(E_H)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
            RCS_HV(index) = real(sum(E_V * conj(E_V)))/real(sum(E_0 * conj(E_0)) * 4.0*M_PI);
          }
          else {
            blitz::Array<std::complex<double>, 1> E_V(3);
            E_V = sum(SAR_local_y_hat * E_obs(0, all));
            RCS_VV(index) = real(sum(E_V * conj(E_V)))/real(sum((E_0 * conj(E_0))) * 4.0*M_PI);
          }
          if (my_id==master) {
            writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HH_ASCII.txt", RCS_HH);
            writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HV_ASCII.txt", RCS_HV);
            writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_VV_ASCII.txt", RCS_VV);
            writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "r_SAR.txt", r_SAR);
          }
        }
      }
    }
  }
  // final writings
  if (my_id==master) {
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HH_ASCII.txt", RCS_HH);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_HV_ASCII.txt", RCS_HV);
    writeFloatBlitzArray1DToASCIIFile(RESULT_DATA_PATH + "SAR_RCS_VV_ASCII.txt", RCS_VV);
    writeFloatBlitzArray2DToASCIIFile(RESULT_DATA_PATH + "r_SAR.txt", r_SAR);
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "numberOfMatvecs.txt", octtree.getNumberOfUpdates());
    writeIntToASCIIFile(ITERATIVE_DATA_PATH + "iter.txt", iter);
    cout << endl;
  }
}

/****************************************************************************/
/******************************* main ***************************************/
/****************************************************************************/

int main(void) {

  MPI::Init();
  int ierror, num_procs = MPI::COMM_WORLD.Get_size(), my_id = MPI::COMM_WORLD.Get_rank();

  // general variables
  const string TMP = "./tmp" + intToString(my_id), OCTTREE_DATA_PATH = TMP + "/octtree_data/", MESH_DATA_PATH = TMP + "/mesh/", V_CFIE_DATA_PATH = TMP + "/V_CFIE/", RESULT_DATA_PATH = "./result/";
  const string ITERATIVE_DATA_PATH = TMP + "/iterative_data/";

  Mesh target_mesh(MESH_DATA_PATH);
  const int N_RWG = target_mesh.E;
  Octtree octtree(OCTTREE_DATA_PATH, target_mesh.cubes_centroids, my_id, num_procs);
  octtree.computeGaussLocatedArguments(target_mesh);
  // every process will receive about the same number of lines
  blitz::Array<int, 1> localIndexesForSolver;
  const int N_RWGs_per_process = static_cast<int>(floor((N_RWG * 1.0)/num_procs));
  int startIndex = my_id * N_RWGs_per_process, stopIndex = startIndex + N_RWGs_per_process - 1;
  if (my_id == num_procs-1) stopIndex = N_RWG - 1;
  const int N_lines = stopIndex - startIndex + 1;
  localIndexesForSolver.resize(N_lines);
  for (int j=0 ; j<N_lines ; ++j) localIndexesForSolver(j) = j + startIndex;
  // we create a local mesh
  LocalMesh local_target_mesh(target_mesh, localIndexesForSolver);
  local_target_mesh.writeLocalMeshToFile(MESH_DATA_PATH);
  target_mesh.~Mesh();
  // OK, what kind of simulation do we want to run?
  octtree.constructArrays();
  int BISTATIC, MONOSTATIC_RCS, MONOSTATIC_SAR;
  readIntFromASCIIFile(TMP + "/BISTATIC.txt", BISTATIC);
  readIntFromASCIIFile(TMP + "/MONOSTATIC_RCS.txt", MONOSTATIC_RCS);
  readIntFromASCIIFile(TMP + "/MONOSTATIC_SAR.txt", MONOSTATIC_SAR);
  string SOLVER;
  readStringFromASCIIFile(ITERATIVE_DATA_PATH + "SOLVER.txt", SOLVER);
  if (my_id==0) cout << "SOLVER IS = " << SOLVER << endl;
  ierror = MPI_Barrier(MPI::COMM_WORLD);
  // bistatic computation
  if (BISTATIC==1) computeForOneExcitation(octtree, local_target_mesh, SOLVER, TMP, OCTTREE_DATA_PATH, MESH_DATA_PATH, V_CFIE_DATA_PATH, RESULT_DATA_PATH, ITERATIVE_DATA_PATH);
  // monostatic RCS computation
  local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
  if (MONOSTATIC_RCS==1) computeMonostaticRCS(octtree, local_target_mesh, SOLVER, TMP, OCTTREE_DATA_PATH, MESH_DATA_PATH, V_CFIE_DATA_PATH, RESULT_DATA_PATH, ITERATIVE_DATA_PATH);
  // monostatic SAR computation
  local_target_mesh.setLocalMeshFromFile(MESH_DATA_PATH);
  if (MONOSTATIC_SAR==1) computeMonostaticSAR(octtree, local_target_mesh, SOLVER, TMP, OCTTREE_DATA_PATH, MESH_DATA_PATH, V_CFIE_DATA_PATH, RESULT_DATA_PATH, ITERATIVE_DATA_PATH);
  MPI::Finalize();
  return 0;
}
