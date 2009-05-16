/***************************************************************************
 * generate_K_G_J_grid.cpp  function that generates 3D tables for the DGFs
 *                          those tables will be interpolated in the MoM
 *
 * Copyright (C) 2000-2005 Idesbald van den Bosch <vandenbosch@emic.ucl.ac.be>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * Suggestions:          vandenbosch@emic.ucl.ac.be
 * Bugs:                 vandenbosch@emic.ucl.ac.be
 *
 * For more information, please see the ... Home Page:
 *
 ****************************************************************************/
#include <iostream>
#include <complex>
#include <blitz/array.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "layers_constants.h"
#include "K_G_J.h"
#include "K_G_grid.h"
#include "splint.h"

void generate_grid_points (Array<double, 1>& rho_values, Array<double, 1>& z_values, Array<double, 1>& z_prime_values, const double rho_min, const double rho_max, const double z_min, const double z_max, const double z_prime_min, const double z_prime_max, const double lambda_layer) {

  firstIndex k;
  double Delta_z = lambda_layer/12, Delta_z_prime = lambda_layer/12, Delta_rho = lambda_layer/12;
  int N_Delta_z = 2 * (int) ceil((z_max-z_min)/(2*Delta_z));
  int N_Delta_z_prime = 2 * (int) ceil((z_prime_max-z_prime_min)/(2*Delta_z_prime));
  int N_Delta_rho = 2 * (int) ceil((rho_max-rho_min)/(2*Delta_rho));
  if ((N_Delta_rho<3) && (rho_min<rho_max)) N_Delta_rho = 3;
  if ((N_Delta_z<3) && (z_min<z_max)) N_Delta_z = 3;
  if ((N_Delta_z_prime<3) && (z_prime_min<z_prime_max)) N_Delta_z_prime = 3;
  int N_z = N_Delta_z+1, N_rho = N_Delta_rho+1, N_z_prime = N_Delta_z_prime+1;
  Delta_z = (z_max-z_min)/N_Delta_z;
  Delta_z_prime = (z_prime_max-z_prime_min)/N_Delta_z_prime;
  Delta_rho = (rho_max-rho_min)/N_Delta_rho;

  // grid points definition
  rho_values.resize(N_rho);
  z_values.resize(N_z);
  z_prime_values.resize(N_z_prime);
  rho_values = rho_min + k*Delta_rho;
  rho_values (N_rho-1) = rho_max;
  z_values = z_min + k*Delta_z;
  z_values (N_z-1) = z_max;
  z_prime_values =  z_prime_min + k*Delta_z_prime;
  z_prime_values (N_z_prime-1) = z_prime_max;
}

void generate_KA (const int USE_SYMMETRY, KA_grid & KA_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC) {

  int m, n, p, N_rho = rho_values.size(), N_z = z_values.size(), N_z_prime = z_prime_values.size(), N_layers = LC.eps_i.size();

  KA_grid KA_tab_coarse;
  KA_tab_coarse.xx.resize(N_z, N_z_prime, N_rho);
  KA_tab_coarse.xz_yz.resize(N_z, N_z_prime, N_rho);
  KA_tab_coarse.zx_zy.resize(N_z, N_z_prime, N_rho);
  KA_tab_coarse.zz.resize(N_z, N_z_prime, N_rho);

  cout << "generation of KA_J : " << N_z << " * " << N_z_prime << " * " << N_rho << " = " << N_z_prime * N_z * N_rho << endl;
  for (p=0 ; p<N_rho ; p++) {
    for (m=0 ; m<N_z ; m++) {
      int ind_z = count(LC.z_i<z_values(m));
      for (n=0 ; n<N_z_prime ; n++) {

        int ind_z_prime = count(LC.z_i<z_prime_values(n));
        if (USE_SYMMETRY==1) { // we use the symmetry in z and z' in the Green's functions
          if (m>=n) { // we compute KA_xx and KA_zz only if (z index >= z' index)
            KA_tab_coarse.xx (m, n, p) = KA_xx (rho_values(p), z_values(m), z_prime_values(n), LC);
            KA_tab_coarse.zz (m, n, p) = KA_zz (rho_values(p), z_values(m), z_prime_values(n), LC);
            KA_tab_coarse.xx (n, m, p) = KA_tab_coarse.xx (m, n, p);
            KA_tab_coarse.zz (n, m, p) = KA_tab_coarse.zz (m, n, p);
	  }
          KA_tab_coarse.zx_zy (m, n, p) = KA_zx_zy (rho_values(p), z_values(m), z_prime_values(n), LC);
          KA_tab_coarse.xz_yz (n, m, p) = -LC.mu_i(ind_z_prime)/LC.mu_i(ind_z) * KA_tab_coarse.zx_zy (m, n, p);
	}
        else { // if (USE_SYMMETRY!=1)
          KA_tab_coarse.xx (m, n, p) = KA_xx (rho_values(p), z_values(m), z_prime_values(n), LC);
          KA_tab_coarse.zz (m, n, p) = KA_zz (rho_values(p), z_values(m), z_prime_values(n), LC);
          KA_tab_coarse.zx_zy (m, n, p) = KA_zx_zy (rho_values(p), z_values(m), z_prime_values(n), LC);
          KA_tab_coarse.xz_yz (m, n, p) = KA_xz_yz (rho_values(p), z_values(m), z_prime_values(n), LC);
	}

      }
    }
  }
  decimate_3D (z_values, z_prime_values, rho_values, KA_tab_coarse.xx, N_decim_z, N_decim_z_prime, N_decim_rho, KA_tab.z_values, KA_tab.z_prime_values, KA_tab.rho_values, KA_tab.xx);
  decimate_3D (z_values, z_prime_values, rho_values, KA_tab_coarse.xz_yz, N_decim_z, N_decim_z_prime, N_decim_rho, KA_tab.z_values, KA_tab.z_prime_values, KA_tab.rho_values, KA_tab.xz_yz);
  decimate_3D (z_values, z_prime_values, rho_values, KA_tab_coarse.zx_zy, N_decim_z, N_decim_z_prime, N_decim_rho, KA_tab.z_values, KA_tab.z_prime_values, KA_tab.rho_values, KA_tab.zx_zy);
  decimate_3D (z_values, z_prime_values, rho_values, KA_tab_coarse.zz, N_decim_z, N_decim_z_prime, N_decim_rho, KA_tab.z_values, KA_tab.z_prime_values, KA_tab.rho_values, KA_tab.zz);
}

void generate_Kphi (const int USE_SYMMETRY, Kphi_grid & Kphi_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC) {

  int m, n, p, N_rho = rho_values.size(), N_z = z_values.size(), N_z_prime = z_prime_values.size(), N_layers = LC.eps_i.size();

  Kphi_grid Kphi_tab_coarse;
  Kphi_tab_coarse.K.resize (N_z, N_z_prime, N_rho);

  cout << "generation of K_phi : " << N_z << " * " << N_z_prime << " * " << N_rho << " = " << N_z_prime * N_z * N_rho << endl;
  for (p=0 ; p<N_rho ; p++) {
    for (m=0 ; m<N_z ; m++) {
      for (n=0 ; n<N_z_prime ; n++) {

        if (USE_SYMMETRY==1) { // we use the symmetry in z and z' in the Green's functions
          if (m>=n) { // we compute K_phi only if (z index >= z' index)
            Kphi_tab_coarse.K (m, n, p) = K_phi (rho_values(p), z_values(m), z_prime_values(n), LC);
            Kphi_tab_coarse.K (n, m, p) = Kphi_tab_coarse.K (m, n, p); 
	  }
	}
        else Kphi_tab_coarse.K (m, n, p) = K_phi (rho_values(p), z_values(m), z_prime_values(n), LC);

      }
    }
  }
  decimate_3D (z_values, z_prime_values, rho_values, Kphi_tab_coarse.K, N_decim_z, N_decim_z_prime, N_decim_rho, Kphi_tab.z_values, Kphi_tab.z_prime_values, Kphi_tab.rho_values, Kphi_tab.K);
}

void generate_grad_Kphi (const int USE_SYMMETRY, Kphi_grid & Kphi_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC) {

  int m, n, p, N_rho = rho_values.size(), N_z = z_values.size(), N_z_prime = z_prime_values.size(), N_layers = LC.eps_i.size();

  Kphi_grid Kphi_tab_coarse;
  Kphi_tab_coarse.grad_x.resize (N_z, N_z_prime, N_rho);
  Kphi_tab_coarse.grad_z.resize (N_z, N_z_prime, N_rho);

  cout << "generation of grad_Kphi : " << N_z << " * " << N_z_prime << " * " << N_rho << " = " << N_z_prime * N_z * N_rho << endl;
  for (p=0 ; p<N_rho ; p++) {
    for (m=0 ; m<N_z ; m++) {
      for (n=0 ; n<N_z_prime ; n++) {

        if (USE_SYMMETRY==1) { // we use the symmetry in z and z' in the Green's functions
          if (m>=n) { // we compute K_phi only if (z index >= z' index)
            Kphi_tab_coarse.grad_x (m, n, p) = grad_K_phi_xy (rho_values(p), z_values(m), z_prime_values(n), LC);
            Kphi_tab_coarse.grad_x (n, m, p) = Kphi_tab_coarse.grad_x (m, n, p);
	  }
	}
	else Kphi_tab_coarse.grad_x (m, n, p) = grad_K_phi_xy (rho_values(p), z_values(m), z_prime_values(n), LC);
        Kphi_tab_coarse.grad_z (m, n, p) = grad_K_phi_z (rho_values(p), z_values(m), z_prime_values(n), LC);
      }
    }
  }
  decimate_3D (z_values, z_prime_values, rho_values, Kphi_tab_coarse.grad_x, N_decim_z, N_decim_z_prime, N_decim_rho, Kphi_tab.z_values, Kphi_tab.z_prime_values, Kphi_tab.rho_values, Kphi_tab.grad_x);
  decimate_3D (z_values, z_prime_values, rho_values, Kphi_tab_coarse.grad_z, N_decim_z, N_decim_z_prime, N_decim_rho, Kphi_tab.z_values, Kphi_tab.z_prime_values, Kphi_tab.rho_values, Kphi_tab.grad_z);
}

void generate_grad_prime_Kphi (const int USE_SYMMETRY, Kphi_grid & Kphi_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC) {

  int m, n, p, N_rho = rho_values.size(), N_z = z_values.size(), N_z_prime = z_prime_values.size(), N_layers = LC.eps_i.size();

  Kphi_grid Kphi_tab_coarse;
  Kphi_tab_coarse.grad_prime_x.resize (N_z, N_z_prime, N_rho);
  Kphi_tab_coarse.grad_prime_z.resize (N_z, N_z_prime, N_rho);

  cout << "generation of grad_prime_Kphi : " << N_z << " * " << N_z_prime << " * " << N_rho << " = " << N_z_prime * N_z * N_rho << endl;
  for (p=0 ; p<N_rho ; p++) {
    for (m=0 ; m<N_z ; m++) {
      for (n=0 ; n<N_z_prime ; n++) {

        if (USE_SYMMETRY==1) { // we use the symmetry in z and z' in the Green's functions
          if (m>=n) { // we compute K_phi only if (z index >= z' index)
            Kphi_tab_coarse.grad_prime_x (m, n, p) = grad_prime_K_phi_xy (rho_values(p), z_values(m), z_prime_values(n), LC);
            Kphi_tab_coarse.grad_prime_x (n, m, p) = Kphi_tab_coarse.grad_prime_x (m, n, p);
	  }
	}
	else Kphi_tab_coarse.grad_prime_x (m, n, p) = grad_prime_K_phi_xy (rho_values(p), z_values(m), z_prime_values(n), LC);
        Kphi_tab_coarse.grad_prime_z (m, n, p) = grad_prime_K_phi_z (rho_values(p), z_values(m), z_prime_values(n), LC);
      }
    }
  }
  decimate_3D (z_values, z_prime_values, rho_values, Kphi_tab_coarse.grad_prime_x, N_decim_z, N_decim_z_prime, N_decim_rho, Kphi_tab.z_values, Kphi_tab.z_prime_values, Kphi_tab.rho_values, Kphi_tab.grad_prime_x);
  decimate_3D (z_values, z_prime_values, rho_values, Kphi_tab_coarse.grad_prime_z, N_decim_z, N_decim_z_prime, N_decim_rho, Kphi_tab.z_values, Kphi_tab.z_prime_values, Kphi_tab.rho_values, Kphi_tab.grad_prime_z);
}

void generate_G_HJ (const int USE_SYMMETRY, G_HJ_grid & G_HJ_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC) {

  int m, n, p, N_rho = rho_values.size(), N_z = z_values.size(), N_z_prime = z_prime_values.size(), N_layers = LC.eps_i.size();

  G_HJ_grid G_HJ_tab_coarse;
  G_HJ_tab_coarse.xx.resize (N_z, N_z_prime, N_rho);
  G_HJ_tab_coarse.xy_1.resize (N_z, N_z_prime, N_rho);
  G_HJ_tab_coarse.xz_yz.resize (N_z, N_z_prime, N_rho);
  G_HJ_tab_coarse.zx_zy.resize (N_z, N_z_prime, N_rho);

  cout << "generation of G_HJ : " << N_z << " * " << N_z_prime << " * " << N_rho << " = " << N_z_prime * N_z * N_rho << endl;
  for (p=0 ; p<N_rho ; p++) {
    for (m=0 ; m<N_z ; m++) {
      for (n=0 ; n<N_z_prime ; n++) {

        G_HJ_tab_coarse.xx (m, n, p) = G_HJ_xx (rho_values(p), z_values(m), z_prime_values(n), LC); // S_2 {I_i_e - I_i_h}
        G_HJ_tab_coarse.xy_1 (m, n, p) = G_HJ_xy_1 (rho_values(p), z_values(m), z_prime_values(n), LC); // S_0 {I_i_h + I_i_e}
        if (USE_SYMMETRY==1) { // we use the symmetry in z and z' in the Green's functions
          if (m>=n) { // we compute G_HJ_xz_yz and G_HJ_zx_zy only if (z index >= z' index)
	    G_HJ_tab_coarse.xz_yz (m, n, p) = G_HJ_xz_yz (rho_values(p), z_values(m), z_prime_values(n), LC); // j/(w*LC.eps_0*LC.eps(n)) * S_1 {k_rho * I_v_e}
            G_HJ_tab_coarse.zx_zy (m, n, p) = G_HJ_zx_zy (rho_values(p), z_values(m), z_prime_values(n), LC); // j/(w*LC.mu_0*LC.mu(m)) * S_1 {k_rho * V_i_h}
            G_HJ_tab_coarse.xz_yz (n, m, p) = G_HJ_tab_coarse.xz_yz (m, n, p);
            G_HJ_tab_coarse.zx_zy (n, m, p) = G_HJ_tab_coarse.zx_zy (m, n, p);
	  }
	}
        else { // if (USE_SYMMETRY!=1)
          G_HJ_tab_coarse.xz_yz (m, n, p) = G_HJ_xz_yz (rho_values(p), z_values(m), z_prime_values(n), LC); // j/(w*LC.eps_0*LC.eps(n)) * S_1 {k_rho * I_v_e}
          G_HJ_tab_coarse.zx_zy (m, n, p) = G_HJ_zx_zy (rho_values(p), z_values(m), z_prime_values(n), LC); // j/(w*LC.mu_0*LC.mu(m)) * S_1 {k_rho * V_i_h}
	}

      }
    }
  }

  decimate_3D (z_values, z_prime_values, rho_values, G_HJ_tab_coarse.xx, N_decim_z, N_decim_z_prime, N_decim_rho, G_HJ_tab.z_values, G_HJ_tab.z_prime_values, G_HJ_tab.rho_values, G_HJ_tab.xx);
  decimate_3D (z_values, z_prime_values, rho_values, G_HJ_tab_coarse.xy_1, N_decim_z, N_decim_z_prime, N_decim_rho, G_HJ_tab.z_values, G_HJ_tab.z_prime_values, G_HJ_tab.rho_values, G_HJ_tab.xy_1);
  decimate_3D (z_values, z_prime_values, rho_values, G_HJ_tab_coarse.xz_yz, N_decim_z, N_decim_z_prime, N_decim_rho, G_HJ_tab.z_values, G_HJ_tab.z_prime_values, G_HJ_tab.rho_values, G_HJ_tab.xz_yz);
  decimate_3D (z_values, z_prime_values, rho_values, G_HJ_tab_coarse.zx_zy, N_decim_z, N_decim_z_prime, N_decim_rho, G_HJ_tab.z_values, G_HJ_tab.z_prime_values, G_HJ_tab.rho_values, G_HJ_tab.zx_zy);
}

void generate_G_EJ (const int USE_SYMMETRY, G_EJ_grid & G_EJ_tab, const Array<double, 1>& z_values, const Array<double, 1>& z_prime_values, const Array<double, 1>& rho_values, const int N_decim_z, const int N_decim_z_prime, const int N_decim_rho, const layers_constants & LC) {

  int m, n, p, N_rho = rho_values.size(), N_z = z_values.size(), N_z_prime = z_prime_values.size(), N_layers = LC.eps_i.size();
  Array<complex<double>, 1> G_EJ_tmp (5);

  G_EJ_grid G_EJ_tab_coarse;
  G_EJ_tab_coarse.xx_1.resize (N_z, N_z_prime, N_rho);
  G_EJ_tab_coarse.xx_2.resize (N_z, N_z_prime, N_rho);
  G_EJ_tab_coarse.xz_yz.resize (N_z, N_z_prime, N_rho);
  G_EJ_tab_coarse.zx_zy.resize (N_z, N_z_prime, N_rho);
  G_EJ_tab_coarse.zz.resize (N_z, N_z_prime, N_rho);

  cout << "generation of G_EJ : " << N_z << " * " << N_z_prime << " * " << N_rho << " = " << N_z_prime * N_z * N_rho << endl;
  for (p=0 ; p<N_rho ; p++) {
    for (m=0 ; m<N_z ; m++) {
      for (n=0 ; n<N_z_prime ; n++) {
        // no use of symmetry as G_EJ (and G_HM) is used only in computation of V
        G_EJ_tmp = G_EJ (rho_values(p), z_values(m), z_prime_values(n), LC);
        G_EJ_tab_coarse.xx_1 (m, n, p) = G_EJ_tmp (0); // S_0 {V_i_e + V_i_h}
        G_EJ_tab_coarse.xx_2 (m, n, p) = G_EJ_tmp (1); // S_2 {V_i_e - V_i_h}
	G_EJ_tab_coarse.xz_yz (m, n, p) = G_EJ_tmp (2); // -j/(w*LC.eps_0*LC.eps(n)) * S_1 {k_rho * V_v_e}
        G_EJ_tab_coarse.zx_zy (m, n, p) = G_EJ_tmp (3); // -j/(w*LC.eps_0*LC.eps(m)) * S_1 {k_rho * I_i_e}
        G_EJ_tab_coarse.zz (m, n, p) = G_EJ_tmp (4); // -1/(w^2*LC.eps_0^2*LC.eps(n)*LC.eps(m)) * S_0 {k_rho^2 * I_v_e}

      }
    }
  }

  decimate_3D (z_values, z_prime_values, rho_values, G_EJ_tab_coarse.xx_1, N_decim_z, N_decim_z_prime, N_decim_rho, G_EJ_tab.z_values, G_EJ_tab.z_prime_values, G_EJ_tab.rho_values, G_EJ_tab.xx_1);
  decimate_3D (z_values, z_prime_values, rho_values, G_EJ_tab_coarse.xx_2, N_decim_z, N_decim_z_prime, N_decim_rho, G_EJ_tab.z_values, G_EJ_tab.z_prime_values, G_EJ_tab.rho_values, G_EJ_tab.xx_2);
  decimate_3D (z_values, z_prime_values, rho_values, G_EJ_tab_coarse.xz_yz, N_decim_z, N_decim_z_prime, N_decim_rho, G_EJ_tab.z_values, G_EJ_tab.z_prime_values, G_EJ_tab.rho_values, G_EJ_tab.xz_yz);
  decimate_3D (z_values, z_prime_values, rho_values, G_EJ_tab_coarse.zx_zy, N_decim_z, N_decim_z_prime, N_decim_rho, G_EJ_tab.z_values, G_EJ_tab.z_prime_values, G_EJ_tab.rho_values, G_EJ_tab.zx_zy);
  decimate_3D (z_values, z_prime_values, rho_values, G_EJ_tab_coarse.zz, N_decim_z, N_decim_z_prime, N_decim_rho, G_EJ_tab.z_values, G_EJ_tab.z_prime_values, G_EJ_tab.rho_values, G_EJ_tab.zz);
}
