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
#include "R_G_soil.h"

void z_ant_min_z_ant_max (double & z_ant_min, double & z_ant_max, const TinyVector<double, 3>& r_ant, const TinyVector<double, 3>& x_hat_ant, const TinyVector<double, 3>& y_hat_ant) {

  int j;

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
  TinyVector<double, 3> r_tmp;
  for (j=0 ; j<nodes_antenna_tmp.rows() ; j++) { // we loop on all the nodes of the antenna mesh
    r_tmp = r_ant + x_hat_ant * nodes_antenna_tmp (j, 1) + y_hat_ant * nodes_antenna_tmp (j, 2);
    z_ant_min = (z_ant_min < r_tmp (2)) ? z_ant_min : r_tmp (2);
    z_ant_max = (z_ant_max > r_tmp (2)) ? z_ant_max : r_tmp (2);
  }
}

void R_G_soil_main (complex<double>& R11_dip_soil, complex<double>& G11_dip_soil, complex<double>& R11_horn_soil, complex<double>& G11_horn_soil) {
  int i, E, USE_SYMMETRY;
  TinyVector<double, 3> r_ant, x_hat_ant, y_hat_ant, r_tmp;
  TinyVector<double, 3> y_hat (0.0, 1.0, 0.0);
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

  TinyVector<double, 3> r_VP;
  r_VP = r_ant;
  r_VP (2) += 0.0714;

  mesh MESH_ANTENNA;
  char * nodes_antenna = "nodes_antenna.txt", * triangles_antenna = "triangles_antenna.txt", * edges_antenna = "edges_antenna.txt";
  create_mesh (MESH_ANTENNA, nodes_antenna, triangles_antenna, edges_antenna);

  // Green function domain definition
  double eps_limit = 1.0e-8;
  double z_ant_min, z_ant_max;
  z_ant_min_z_ant_max (z_ant_min, z_ant_max, r_ant, x_hat_ant, y_hat_ant);

  layers_constants LC;
  create_LC (LC, eps_0, mu_0);
  int ind_layer;
  if (LC.N>1) {
    ind_layer = count(LC.z_i<z_ant_min);
    if (ind_layer!=count(LC.z_i<z_ant_max)) cout << "BIG PROBLEM : ANTENNA CROSSES INTERFACES!!" << endl;
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

  double lambda_0 = c/LC.w * 2.0 * M_PI;
  double J_0 = sqrt(lambda_0/sqrt(mu_0/eps_0) * 3.0/M_PI); // for P_rad = 1 \cite[eq. 4-16]{Balanis_82};
  J_ant = -J_0 * y_hat_ant;
  M_ant = 0.0 * x_hat_ant;

  Array<double, 1> rho_values, z_values, z_prime_values;

  // generation of dipolar soil response.
  cout << endl << "/********************************************************************************/" << endl;
  cout <<         "/*         computation of soil response due to DIPOLE excitation                */" << endl;
  cout <<         "/********************************************************************************/" << endl;

  generate_grid_points (rho_values, z_values, z_prime_values, MESH_ANTENNA.rho_min, MESH_ANTENNA.rho_max/3.0, r_VP (2), r_VP (2), r_VP (2), r_VP (2), lambda_layer);
  cout << "antenna rho max/3.0 = " << MESH_ANTENNA.rho_max/3.0 << endl;

  USE_SYMMETRY = 1;
  int N_decim_z = (z_values.size()>3) ? 1 : 0, N_decim_z_prime = (z_prime_values.size()>3) ? 1 : 0, N_decim_rho = (rho_values.size()>3) ? 1 : 0;
  cout << "Generation of multilayered media green functions (antenna to antenna) for J currents on dipole" << endl;

  G_EJ_grid G_EJ_tab_VP_VP;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_VP_VP, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);
  G_HJ_grid G_HJ_tab_VP_VP;
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_VP_VP, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);
  TinyVector<complex<double>, 3> E_dip, H_dip;
  E_H_J_dip (E_dip, H_dip, G_EJ_tab_VP_VP, G_HJ_tab_VP_VP, J_ant, r_VP, r_VP, LC);

  G11_dip_soil = dot (E_dip, y_hat); // the field is projected on y_hat
  R11_dip_soil = -0.5 * dot (E_dip, J_ant); // -1/2 e_{soil, 1} * j_1

  // generation of cosine soil response.
  cout << endl << "/********************************************************************************/" << endl;
  cout <<         "/*         computation of soil R11 due to COSINE excitation                     */" << endl;
  cout <<         "/********************************************************************************/" << endl;
  z_ant_min_z_ant_max (z_ant_min, z_ant_max, r_ant, x_hat_ant, y_hat_ant);

  generate_grid_points (rho_values, z_values, z_prime_values, MESH_ANTENNA.rho_min, MESH_ANTENNA.rho_max, z_ant_min, z_ant_max, z_ant_min, z_ant_max, lambda_layer);
  cout << "antenna rho max = " << MESH_ANTENNA.rho_max << ", antenna z min = " << z_ant_min << ", antenna z max = " << z_ant_max << endl;

  USE_SYMMETRY = 1;
  N_decim_z = (z_values.size()>3) ? 1 : 0, N_decim_z_prime = (z_prime_values.size()>3) ? 1 : 0, N_decim_rho = (rho_values.size()>3) ? 1 : 0;

  cout << "Generation of multilayered media green functions (antenna to antenna) for J currents on aperture" << endl;
  G_EJ_grid G_EJ_tab_ant_ant;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_ant_ant, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);
  G_HJ_grid G_HJ_tab_ant_ant;
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_ant_ant, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);

  cout << "Generation of multilayered media dual green functions (antenna to antenna) for M currents on antenna" << endl;
  G_EJ_grid G_EJ_tab_ant_ant_dual;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_ant_ant_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);
  G_HJ_grid G_HJ_tab_ant_ant_dual;
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_ant_ant_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);

  R11_horn_soil = R_horn_soil (MESH_ANTENNA, r_ant, x_hat_ant, y_hat_ant, G_EJ_tab_ant_ant, G_HJ_tab_ant_ant, LC, G_EJ_tab_ant_ant_dual, G_HJ_tab_ant_ant_dual, LC_dual);

  cout << endl << "/********************************************************************************/" << endl;
  cout <<         "/*         computation of soil G11 due to COSINE excitation                     */" << endl;
  cout <<         "/********************************************************************************/" << endl;
  z_ant_min_z_ant_max (z_ant_min, z_ant_max, r_ant, x_hat_ant, y_hat_ant);

  generate_grid_points (rho_values, z_values, z_prime_values, MESH_ANTENNA.rho_min, MESH_ANTENNA.rho_max, r_VP(2), r_VP(2), z_ant_min, z_ant_max, lambda_layer);
  USE_SYMMETRY = 0;
  N_decim_z = (z_values.size()>3) ? 1 : 0, N_decim_z_prime = (z_prime_values.size()>3) ? 1 : 0, N_decim_rho = (rho_values.size()>3) ? 1 : 0;

  cout << "Generation of multilayered media green functions (antenna to antenna) for J currents on aperture" << endl;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_ant_ant, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_ant_ant, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC);

  cout << "Generation of multilayered media dual green functions (antenna to antenna) for M currents on antenna" << endl;
  generate_G_EJ (USE_SYMMETRY, G_EJ_tab_ant_ant_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);
  generate_G_HJ (USE_SYMMETRY, G_HJ_tab_ant_ant_dual, z_values, z_prime_values, rho_values, N_decim_z, N_decim_z_prime, N_decim_rho, LC_dual);

  TinyVector<complex<double>, 3> E_J_M, H_J_M;
  E_H_J_M_cosine (E_J_M, H_J_M, MESH_ANTENNA, r_ant, x_hat_ant, y_hat_ant, r_VP, G_EJ_tab_ant_ant, G_HJ_tab_ant_ant, LC, G_EJ_tab_ant_ant_dual, G_HJ_tab_ant_ant_dual, LC_dual);
  G11_horn_soil = dot(E_J_M, y_hat);
  
}
