#ifndef V_E_V_H_H
#define V_E_V_H_H

using namespace std;

#include "mesh.h"

void G_EJ_G_HJ (blitz::Array<std::complex<double>, 2>& G_EJ, 
                blitz::Array<std::complex<double>, 2>& G_HJ, 
                const double r_dip[], 
                const double r_obs[], 
                const std::complex<double>& eps, 
                const std::complex<double>& mu, 
                const std::complex<double>& k);

void V_EJ_HJ_dipole (blitz::Array<std::complex<double>, 1> V_tE_J,
                     blitz::Array<std::complex<double>, 1> V_nE_J,
                     blitz::Array<std::complex<double>, 1> V_tH_J,
                     blitz::Array<std::complex<double>, 1> V_nH_J,
                     const blitz::Array<std::complex<double>, 1>& J_dip,
                     const blitz::Array<double, 1>& r_dip,
                     const blitz::Array<int, 1>& numbers_RWG_test,
                     const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                     const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                     const blitz::Array<double, 2>& RWGNumber_vertexesCoord,
                     const blitz::Array<double, 2>& RWGNumber_oppVertexesCoord,
                     const double w,
                     const std::complex<double>& eps_r,
                     const std::complex<double>& mu_r,
                     const int FULL_PRECISION);

void V_CFIE_dipole (blitz::Array<std::complex<float>, 1> V_CFIE,
                    const blitz::Array<std::complex<float>, 1>& CFIE,
                    const blitz::Array<std::complex<double>, 1>& J_dip,
                    const blitz::Array<double, 1>& r_dip,
                    const blitz::Array<int, 1>& numbers_RWG_test,
                    const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                    const blitz::Array<float, 2>& RWGNumber_trianglesCoord,
                    const double w,
                    const std::complex<double>& eps_r,
                    const std::complex<double>& mu_r,
                    const int FULL_PRECISION);

void local_V_CFIE_dipole (blitz::Array<std::complex<float>, 1>& V_CFIE,
                          const blitz::Array<std::complex<double>, 1>& J_dip,
                          const blitz::Array<double, 1>& r_dip,
                          const LocalMesh & local_target_mesh,
                          const double w,
                          const std::complex<double>& eps_r,
                          const std::complex<double>& mu_r,
                          const blitz::Array<std::complex<float>, 1>& CFIE,
                          const int FULL_PRECISION);

void local_V_CFIE_dipole_array (blitz::Array<std::complex<float>, 1>& V_CFIE,
                                const blitz::Array<std::complex<double>, 2>& J_dip,
                                const blitz::Array<double, 2>& r_dip,
                                const LocalMesh & local_target_mesh,
                                const double w,
                                const std::complex<double>& eps_r,
                                const std::complex<double>& mu_r,
                                const blitz::Array<std::complex<float>, 1>& CFIE,
                                const char CURRENT_TYPE,
                                const int FULL_PRECISION);

void V_EJ_HJ_dipole_alternative (blitz::Array<std::complex<double>, 1> V_tE_J, 
                                 blitz::Array<std::complex<double>, 1> V_nE_J,
                                 blitz::Array<std::complex<double>, 1> V_tH_J,
                                 blitz::Array<std::complex<double>, 1> V_nH_J,
                                 const blitz::Array<std::complex<double>, 1>& J_dip,
                                 const blitz::Array<double, 1>& r_dip,
                                 const blitz::Array<int, 1>& numbers_RWG_test,
                                 const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                                 const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                                 const blitz::Array<double, 2>& RWGNumber_vertexesCoord,
                                 const blitz::Array<double, 2>& RWGNumber_oppVertexesCoord,
                                 const double w,
                                 const std::complex<double>& eps_r,
                                 const std::complex<double>& mu_r,
                                 const int FULL_PRECISION);

void V_EJ_HJ_plane (blitz::Array<std::complex<double>, 1> V_tE_J,
                    blitz::Array<std::complex<double>, 1> V_nE_J,
                    blitz::Array<std::complex<double>, 1> V_tH_J,
                    blitz::Array<std::complex<double>, 1> V_nH_J,
                    const blitz::Array<std::complex<double>, 1>& E_0,
                    const blitz::Array<double, 1>& k_hat,
                    const blitz::Array<double, 1>& r_ref,
                    const blitz::Array<int, 1>& numbers_RWG_test,
                    const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                    const blitz::Array<int, 2>& RWGNumber_signedTriangles,
                    const blitz::Array<double, 2>& RWGNumber_vertexesCoord,
                    const blitz::Array<double, 2>& RWGNumber_oppVertexesCoord,
                    const double w,
                    const std::complex<double>& eps_r,
                    const std::complex<double>& mu_r,
                    const int FULL_PRECISION);

void V_CFIE_plane (blitz::Array<std::complex<float>, 1> V_CFIE,
                   const blitz::Array<std::complex<float>, 1>& CFIE,
                   const blitz::Array<std::complex<double>, 1>& E_0,
                   const blitz::Array<double, 1>& k_hat,
                   const blitz::Array<double, 1>& r_ref,
                   const blitz::Array<int, 1>& numbers_RWG_test,
                   const blitz::Array<int, 1>& RWGNumber_CFIE_OK,
                   const blitz::Array<float, 2>& RWGNumber_trianglesCoord,
                   const double w,
                   const std::complex<double>& eps_r,
                   const std::complex<double>& mu_r,
                   const int FULL_PRECISION);

void local_V_CFIE_plane (blitz::Array<std::complex<float>, 1>& V_CFIE,
                         const blitz::Array<std::complex<double>, 1>& E_0,
                         const blitz::Array<double, 1>& k_hat,
                         const blitz::Array<double, 1>& r_ref,
                         const LocalMesh & local_target_mesh,
                         const double w,
                         const std::complex<double>& eps_r,
                         const std::complex<double>& mu_r,
                         const blitz::Array<std::complex<float>, 1>& CFIE,
                         const int FULL_PRECISION);

void local_V_CFIE_slot (blitz::Array<std::complex<float>, 1>& V_CFIE,
                        const std::complex<double> E_0,
                        const blitz::Array<double, 1>& l_hat,
                        const blitz::Array<double, 1>& r_ref,
                        const double slot_length,
                        const LocalMesh & local_target_mesh,
                        const double w,
                        const std::complex<double>& eps_r,
                        const std::complex<double>& mu_r,
                        const blitz::Array<std::complex<float>, 1>& CFIE,
                        const int FULL_PRECISION);

//void V_E_V_H_horn_BBHA (Array<complex<double>,2>& V_EJ, Array<complex<double>,2>& V_HJ, Array<complex<double>,2>& V_EM, Array<complex<double>,2>& V_HM, const mesh & MESH_TARGET, const mesh & MESH_ANTENNA, const Vector<double,3>& r_ant, const Vector<double,3>& x_hat_ant, const Vector<double,3>& y_hat_ant, const G_EJ_grid & G_EJ_tab_ant_object, const G_HJ_grid & G_HJ_tab_ant_object, const layers_constants & LC, const G_EJ_grid & G_EJ_tab_ant_object_dual, const G_HJ_grid & G_HJ_tab_ant_object_dual, const layers_constants & LC_dual);

#endif
