#include <iostream>
#include <fstream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);
const double mu_0 = 4*M_PI*1e-7;
const double c = 2.99792458e8;
const double eps_0 = 1.0/(c*c*mu_0);

#include "layers_constants.h"
#include "create_LC.h"
#include "K_G_grid.h"
#include "mesh.h"
#include "create_mesh.h"
#include "generate_K_G_J_grid.h"
#include "interp_K_G_grid.h"
#include "V_E_V_H.h"

void rho_min_rho_max (double & rho_min, double & rho_max, double & z_ant_min, double & z_ant_max, const TinyVector<double, 3>& r_ant, const TinyVector<double, 3>& x_hat_ant, const TinyVector<double, 3>& y_hat_ant) {

  int j;
  // reading nodes_target.txt into nodes_target_tmp
  Array<double, 2> nodes_target_tmp;
  {
    char * filename = "nodes_target.txt";
    ifstream ifs(filename);
    if (! ifs.is_open()) { 
      cout << "Vmain : Error opening file : " << filename << endl; 
      exit (1); 
    }
    ifs >> nodes_target_tmp;
    ifs.close();
  }
  double x_min = min (nodes_target_tmp (Range::all(), 1)), x_max = max (nodes_target_tmp (Range::all(), 1));
  double y_min = min (nodes_target_tmp (Range::all(), 2)), y_max = max (nodes_target_tmp (Range::all(), 2));

  // reading nodes_antenna.txt into nodes_antenna_tmp
  Array<double, 2> nodes_antenna_tmp;
  {
    char * filename = "nodes_antenna.txt";
    ifstream ifs(filename);
    if (! ifs.is_open()) { 
      cout << "Vmain : Error opening file : " << filename << endl; 
      exit (1); 
    }
    ifs >> nodes_antenna_tmp;
    ifs.close();
  }

  z_ant_min = r_ant (2);
  z_ant_max = r_ant (2);
  rho_max = max (sqrt (pow2 (nodes_target_tmp (Range::all(), 1)-r_ant (0)) + pow2 (nodes_target_tmp (Range::all(), 2)-r_ant (1)))); 
  rho_min = min (sqrt (pow2 (nodes_target_tmp (Range::all(), 1)-r_ant (0)) + pow2 (nodes_target_tmp (Range::all(), 2)-r_ant (1))));
  if ( ((r_ant (0)>=x_min) && (r_ant (0)<=x_max)) && ((r_ant (1)>=y_min) && (r_ant (1)<=y_max)) ) rho_min = 0.0;
  TinyVector<double, 3> r_tmp;
  for (j=0 ; j<nodes_antenna_tmp.rows() ; j++) { // we loop on all the nodes of the antenna mesh
    r_tmp = r_ant + x_hat_ant * nodes_antenna_tmp (j, 1) + y_hat_ant * nodes_antenna_tmp (j, 2);
    z_ant_min = (z_ant_min < r_tmp (2)) ? z_ant_min : r_tmp (2);
    z_ant_max = (z_ant_max > r_tmp (2)) ? z_ant_max : r_tmp (2);
    double rho_max_tmp = max (sqrt (pow2 (nodes_target_tmp (Range::all(), 1)-r_tmp (0)) + pow2 (nodes_target_tmp (Range::all(), 2)-r_tmp (1)))); 
    double rho_min_tmp = min (sqrt (pow2 (nodes_target_tmp (Range::all(), 1)-r_tmp (0)) + pow2 (nodes_target_tmp (Range::all(), 2)-r_tmp (1))));
    rho_max = (rho_max>rho_max_tmp) ? rho_max : rho_max_tmp;
    rho_min = (rho_min<rho_min_tmp) ? rho_min : rho_min_tmp;
    if ( ((r_tmp (0)>=x_min) && (r_tmp (0)<=x_max)) && ((r_tmp (1)>=y_min) && (r_tmp (1)<=y_max)) ) rho_min = 0.0;
  }
  rho_max = rho_max + 0.05*rho_max;
  rho_min = (rho_min>0.0) ? rho_min-0.05*rho_min : 0.0;
  if (rho_min<0.0) rho_min = 0.0;
}

void Vmain (Array<complex<double>, 2>& V_EJ, Array<complex<double>, 2>& V_HJ, Array<complex<double>, 2>& V_EM, Array<complex<double>, 2>& V_HM, Array<complex<double>, 2>& V_EJ_dip, Array<complex<double>, 2>& V_HJ_dip, Array<complex<double>, 2>& V_EM_dip, Array<complex<double>, 2>& V_HM_dip) {
  int i, E;
  TinyVector<double, 3> r_ant, x_hat_ant, y_hat_ant, r_tmp;
  Array<double, 2> ant_ref_frame (3, 3);
  TinyVector<complex<double>, 3> J_ant, M_ant; 

  // reading ant_ref_frame.txt into r_ant, x_hat_ant, y_hat_ant
  {
    char * filename = "ant_ref_frame.txt";
    ifstream ifs(filename);
    if (! ifs.is_open()) { 
      cout << "Vmain : Error opening file : " << filename << endl; 
      exit (1); 
    }
    else ifs >> ant_ref_frame;
    ifs.close();
    for (i=0 ; i<3 ; i++) {
      r_ant (i) = ant_ref_frame (0, i);
      x_hat_ant (i) = ant_ref_frame (1, i);
      y_hat_ant (i) = ant_ref_frame (2, i);
    }
  }

  mesh MESH_TARGET;
  char * nodes_target = "nodes_target.txt", * triangles_target = "triangles_target.txt", * edges_target = "edges_target.txt";
  create_mesh (MESH_TARGET, nodes_target, triangles_target, edges_target);
  
  mesh MESH_ANTENNA;
  char * nodes_antenna = "nodes_antenna.txt", * triangles_antenna = "triangles_antenna.txt", * edges_antenna = "edges_antenna.txt";
  create_mesh (MESH_ANTENNA, nodes_antenna, triangles_antenna, edges_antenna);

  // Green function domain definition
  double eps_limit = 1.0e-8;
  double z_min = MESH_TARGET.z_min-eps_limit, z_max = MESH_TARGET.z_max+eps_limit;
  double z_ant_min, z_ant_max;
  double rho_max, rho_min;
  rho_min_rho_max (rho_min, rho_max, z_ant_min, z_ant_max, r_ant, x_hat_ant, y_hat_ant);
  rho_max += eps_limit;
  if (rho_min !=0.0) rho_min -= eps_limit;
  if (rho_min < 0.0) rho_min = 0.0;

  layers_constants LC;
  create_LC (LC, eps_0, mu_0);
  int ind_layer;
  if (LC.N>1) {
    ind_layer = count(LC.z_i<z_min);
    if (ind_layer!=count(LC.z_i<z_max)) cout << "BIG PROBLEM : OBJECT CROSSES INTERFACES!!" << endl;
  }
  else {
    cout << "Only one layer of medium: impossible to compute V" << endl; 
    exit (1);
  } 
  double lambda_layer = 2.0*M_PI/real(LC.k_i(ind_layer));

  layers_constants LC_dual;
  create_LC (LC_dual, mu_0, eps_0);
  LC_dual.eps_i = LC.mu_i;
  LC_dual.mu_i = LC.eps_i;

  Array<double, 1> rho_values, z_values, z_prime_values;
  generate_grid_points (rho_values, z_values, z_prime_values, rho_min, rho_max, z_min, z_max, z_ant_min, z_ant_max, lambda_layer);
  // cout << rho_values << endl;

  // generation of horn antenna excitation vector
  cout << endl << "/********************************************************************************/" << endl;
  cout <<         "/*       computation of target excitation vector due to COSINE excitation       */" << endl;
  cout <<         "/********************************************************************************/" << endl;
  cout << "rho_min = " << rho_min << ", rho_max = " << rho_max << ", z_ant_min = " << z_ant_min << ", z_ant_max = " << z_ant_max << endl;
  cout << "r_ant = " << r_ant << ", x_hat_ant = " << x_hat_ant << ", y_hat_ant = " << y_hat_ant << endl;
  int USE_SYMMETRY = 0;
  int N_decim_z = (z_values.size()>3) ? 1 : 0, N_decim_z_prime = (z_prime_values.size()>3) ? 1 : 0, N_decim_rho = (rho_values.size()>3) ? 1 : 0;
  cout << "Generation of multilayered media green functions (antenna to target) for J currents on antenna" << endl;

  G_EJ_grid G_EJ_tab_ant_object;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_ant_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);

  G_HJ_grid G_HJ_tab_ant_object;
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_ant_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);

  cout << "Generation of multilayered media dual green functions (antenna to target) for M currents on antenna" << endl;
  G_EJ_grid G_EJ_tab_ant_object_dual;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_ant_object_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);

  G_HJ_grid G_HJ_tab_ant_object_dual;
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_ant_object_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);

  E = MESH_TARGET.edges.rows()/2;
  V_EJ.resize (E, 3); // [f_m   n_X_f_m   comparison_with_other_method] = 1 line of V matrix
  V_HJ.resize (E, 3);
  V_EM.resize (E, 3);
  V_HM.resize (E, 3);
  V_E_V_H_horn_BBHA (V_EJ, V_HJ, V_EM, V_HM, MESH_TARGET, MESH_ANTENNA, r_ant, x_hat_ant, y_hat_ant, G_EJ_tab_ant_object, G_HJ_tab_ant_object, LC, G_EJ_tab_ant_object_dual, G_HJ_tab_ant_object_dual, LC_dual);

  // comparison with target dipole excitation
  cout << endl << "/********************************************************************************/" << endl;
  cout <<         "/*       computation of target excitation vector due to DIPOLE excitation       */" << endl;
  cout <<         "/********************************************************************************/" << endl;
  TinyVector<double, 3> r_VP;
  r_VP = r_ant;
  r_VP (2) += 0.0714;
  V_EJ_dip.resize (E, 3);
  V_HJ_dip.resize (E, 3);
  V_EM_dip.resize (E, 3);
  V_HM_dip.resize (E, 3);

  rho_min_rho_max (rho_min, rho_max, z_ant_min, z_ant_max, r_VP, x_hat_ant, y_hat_ant);
  generate_grid_points (rho_values, z_values, z_prime_values, rho_min, rho_max, z_min, z_max, r_VP (2), r_VP (2), lambda_layer);
  N_decim_z = (z_values.size()>3) ? 1 : 0, N_decim_z_prime = (z_prime_values.size()>3) ? 1 : 0, N_decim_rho = (rho_values.size()>3) ? 1 : 0;
  cout << "Generation of multilayered media green functions (antenna to target) for J currents on dipole" << endl;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_ant_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_ant_object, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);
  double lambda_0 = c/LC.w * 2.0 * M_PI;
  double J_0 = sqrt(lambda_0/sqrt(mu_0/eps_0) * 3.0/M_PI); // for P_rad = 1 \cite[eq. 4-16]{Balanis_82};
  J_ant = -J_0 * y_hat_ant;
  M_ant = 0.0 * x_hat_ant;

  V_dipole_body (V_EJ_dip, V_HJ_dip, MESH_TARGET, G_EJ_tab_ant_object, G_HJ_tab_ant_object, J_ant, r_VP, LC);
  // V_dipole_body (V_HM_dip, V_EM_dip, MESH_TARGET, G_EJ_tab_ant_object_dual, G_HJ_tab_ant_object_dual, M_ant, r_VP, LC_dual);
  V_HM_dip = 0.0;
  V_EM_dip = 0.0;
  
}
