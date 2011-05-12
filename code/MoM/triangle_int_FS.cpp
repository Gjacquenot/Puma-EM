#include <iostream>
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);

#include "GK_triangle.h"
#include "GL.h"
#include "triangle_int.h"
#include "mesh.h"

/****************************************************************************/
/********************************* Triangle *********************************/
/****************************************************************************/

Triangle::Triangle(void)
{
  r_nodes.resize (3);
  m_i_hat.resize (3);
  s_i_hat.resize (3);
}

Triangle::Triangle(const blitz::TinyVector<double, 3>& r0,
                   const blitz::TinyVector<double, 3>& r1,
                   const blitz::TinyVector<double, 3>& r2,
                   const int tr_number)
{
  number = tr_number;
  r_nodes.resize (3);
  m_i_hat.resize (3);
  s_i_hat.resize (3);
  r_nodes(0) = r0;
  r_nodes(1) = r1;
  r_nodes(2) = r2;
  // gravity center
  TinyVector<double, 3> r2_r1, r3_r1, r3_r2, r_grav_r1, r_grav_r2, r_grav_r3;
  r_grav = (r_nodes(0) + r_nodes(1) + r_nodes(2))/3.0;
  r2_r1 = r_nodes(1) - r_nodes(0);
  r3_r2 = r_nodes(2) - r_nodes(1);
  r3_r1 = r_nodes(2) - r_nodes(0);

  // T.n_hat and T.A construction
  n_hat = cross(r2_r1, r3_r1);
  A = sqrt(dot(n_hat, n_hat))/2.0;
  n_hat = n_hat/(2.0*A);

  // T.s_i_hat construction
  s_i_hat (0) = r2_r1/sqrt (dot (r2_r1, r2_r1));
  s_i_hat (1) = r3_r2/sqrt (dot (r3_r2, r3_r2));
  s_i_hat (2) = -1.0*r3_r1/sqrt (dot (r3_r1, r3_r1));

  // T.m_i_hat construction
  for (int i=0 ; i<3 ; i++) m_i_hat(i) = cross(s_i_hat(i), n_hat);

  // R_max computation. R_max is the greatest |r_grav - ri|, i=0..2  
  r_grav_r1 = r_grav - r_nodes(0);
  r_grav_r2 = r_grav - r_nodes(1);
  r_grav_r3 = r_grav - r_nodes(2);
  double R1 = sqrt( r_grav_r1(0) * r_grav_r1(0) + r_grav_r1(1) * r_grav_r1(1) + r_grav_r1(2) * r_grav_r1(2) );
  double R2 = sqrt( r_grav_r2(0) * r_grav_r2(0) + r_grav_r2(1) * r_grav_r2(1) + r_grav_r2(2) * r_grav_r2(2) );
  double R3 = sqrt( r_grav_r3(0) * r_grav_r3(0) + r_grav_r3(1) * r_grav_r3(1) + r_grav_r3(2) * r_grav_r3(2) );
  R_max = ( R1 > R2 ) ? R1 : R2;
  R_max = ( R_max > R3) ? R_max : R3;
}

void Triangle::copyTriangle (const Triangle& triangleToCopy) // copy member function
{
  number = triangleToCopy.number;
  r_nodes.resize (3);
  m_i_hat.resize (3);
  s_i_hat.resize (3);
  for (int i=0 ; i<3 ; ++i) {
    r_nodes(i) = triangleToCopy.r_nodes(i);
    m_i_hat(i) = triangleToCopy.m_i_hat(i);
    s_i_hat(i) = triangleToCopy.s_i_hat(i);
  }
  n_hat = triangleToCopy.n_hat;
  r_grav = triangleToCopy.r_grav;
  A = triangleToCopy.A;
  R_max = triangleToCopy.R_max;
  RWGIndexes.resize(triangleToCopy.RWGIndexes.size());
  RWGIndexes = triangleToCopy.RWGIndexes;
  signInRWG.resize(triangleToCopy.signInRWG.size());
  signInRWG = triangleToCopy.signInRWG;
  indexesInRWGs.resize(triangleToCopy.indexesInRWGs.size());
  indexesInRWGs = triangleToCopy.indexesInRWGs;
}

Triangle::Triangle(const Triangle& triangleToCopy) // copy constructor
{
  copyTriangle(triangleToCopy);
}

Triangle& Triangle::operator=(const Triangle& triangleToCopy) { // copy assignment
  copyTriangle(triangleToCopy);
  return *this;
}

Triangle::~Triangle() {
  RWGIndexes.clear();
  indexesInRWGs.clear();
  signInRWG.clear();
  r_nodes.free();
  m_i_hat.free();
  s_i_hat.free();
}

//! this function is useful in the MoM
void constructVectorTriangles(std::vector<Triangle>& triangles,
                              const std::vector<RWG>& vectorRWGs,
                              const std::vector<Dictionary<int, int> >& TriangleToRWG)
{
  blitz::Range all = blitz::Range::all();
  int index = 0;
  for (int i=0 ; i<TriangleToRWG.size() ; i++) {
    const int tr_number = TriangleToRWG[i].getKey();
    const int RWG_index = TriangleToRWG[i].getVal();
    int indexInRWG = (tr_number == vectorRWGs[RWG_index].triangleNumbers(0)) ? 0 : 1;
    double sign = static_cast<double>(vectorRWGs[RWG_index].triangleSigns(indexInRWG));
    if (triangles.size()==0) { // initialisation of the vector of triangles
      if (indexInRWG==0) triangles.push_back(Triangle(vectorRWGs[RWG_index].vertexesCoord(0),
                                                      vectorRWGs[RWG_index].vertexesCoord(1),
                                                      vectorRWGs[RWG_index].vertexesCoord(2),
                                                      tr_number));
      else triangles.push_back(Triangle(vectorRWGs[RWG_index].vertexesCoord(2),
                                        vectorRWGs[RWG_index].vertexesCoord(1),
                                        vectorRWGs[RWG_index].vertexesCoord(3),
                                        tr_number));
    }
    else if (tr_number!=triangles[index].number) { // we create a new triangle
      index++;
      if (indexInRWG==0) triangles.push_back(Triangle(vectorRWGs[RWG_index].vertexesCoord(0),
                                                      vectorRWGs[RWG_index].vertexesCoord(1),
                                                      vectorRWGs[RWG_index].vertexesCoord(2),
                                                      tr_number));
      else triangles.push_back(Triangle(vectorRWGs[RWG_index].vertexesCoord(2),
                                        vectorRWGs[RWG_index].vertexesCoord(1),
                                        vectorRWGs[RWG_index].vertexesCoord(3),
                                        tr_number));
    }
    // we have to actualise RWGIndexes, signInRWG and indexesInRWGs
    triangles[index].RWGIndexes.push_back(RWG_index);
    triangles[index].signInRWG.push_back(sign);
    triangles[index].indexesInRWGs.push_back(indexInRWG);
  }
  std::vector<Triangle>(triangles).swap(triangles);
}

/****************************************************************************/
/****************************** RWG_function ********************************/
/****************************************************************************/

//! Half_RWG_function class
RWG::RWG(const int RWG_number,
         const blitz::TinyVector<int, 2>& triangle_numbers,
         const blitz::TinyVector<int, 2>& triangle_signs,
         const blitz::Array<double, 1>& r0,
         const blitz::Array<double, 1>& r1,
         const blitz::Array<double, 1>& r2,
         const blitz::Array<double, 1>& r3)
{
  number = RWG_number;
  triangleNumbers = triangle_numbers;
  triangleSigns = triangle_signs;
  vertexesCoord.resize(4);
  for (int j=0 ; j<3 ; j++) {
    vertexesCoord(0)(j) = r0(j);
    vertexesCoord(1)(j) = r1(j);
    vertexesCoord(2)(j) = r2(j);
    vertexesCoord(3)(j) = r3(j);
  }
  blitz::TinyVector<double, 3> r1_r2(vertexesCoord(1)-vertexesCoord(2));
  length = std::sqrt(blitz::dot(r1_r2, r1_r2));
}

void RWG::copyRWG(const RWG& RWGToCopy)
{
  number = RWGToCopy.number;
  triangleNumbers = RWGToCopy.triangleNumbers;
  triangleSigns = RWGToCopy.triangleSigns;
  length = RWGToCopy.length;
  vertexesCoord.resize(4);
  for (int i=0 ; i<4 ; i++) vertexesCoord(i) = RWGToCopy.vertexesCoord(i);
}

RWG::RWG(const RWG& RWGToCopy) // copy constructor
{
  copyRWG(RWGToCopy);
}

RWG& RWG::operator=(const RWG& RWGToCopy) { // copy assignment
  copyRWG(RWGToCopy);
  return *this;
}

/****************************************************************************/
/***************** various integration functions ****************************/
/****************************************************************************/

void IT_fm_fn (double & IT_r_square,
               blitz::TinyVector<double, 3>& IT_r,
               const Triangle & T)
{
  IT_r = T.A * (T.r_nodes (0) + T.r_nodes (1) + T.r_nodes (2))/3.0;
  IT_r_square = T.A * (dot(T.r_nodes (0), T.r_nodes (0)) + dot(T.r_nodes (1), T.r_nodes (1)) + dot(T.r_nodes (2), T.r_nodes (2)) + dot(T.r_nodes (0), T.r_nodes (1)) + dot(T.r_nodes (0), T.r_nodes (2)) + dot(T.r_nodes (1), T.r_nodes (2)))/6.0;
}

void IT_singularities (double & IT_1_R,
                       double & IT_R,
                       blitz::TinyVector<double, 3>& IT_1_R_rprime_r,
                       blitz::TinyVector<double, 3>& IT_R_rprime_r,
                       blitz::TinyVector<double, 3>& IT_grad_1_R,
                       blitz::TinyVector<double, 3>& IT_grad_R,
                       const blitz::TinyVector<double,3>& r,
                       const Triangle & T)
/**
 * notations are taken from I.Hanninen, M. Taskinen and J. Sarvas, "Singularity subtraction integral Formulae
 * for surface integral equations with RWG, rooftop and hybrid basis functions", PIER 63, 243--278, 2006
 */
{
  const double h = dot(r-T.r_nodes (0), T.n_hat), abs_h = abs(h), sign_h = (abs_h<1.0e-10) ? 0.0 : h/abs_h;
  double t_i_0, s_plus__i, s_minus__i, R_plus__i, R_minus__i, R_i_0_square;
  double I_L_minus_1__i, I_L_plus_1__i, I_L_plus_3__i, beta_i, K_1_minus_1__i, K_1_plus_1__i;
  const blitz::TinyVector<double, 3> rho(r - T.n_hat * dot(r, T.n_hat));
  blitz::TinyVector<double, 3> r_plus__i, r_minus__i;
  blitz::TinyVector<double, 3> K_2_minus_1__i, K_2_plus_1__i, K_3_minus_1__i, K_3_plus_1__i;

  IT_1_R = 0.0;
  IT_R = 0.0;
  IT_1_R_rprime_r = 0.0;
  IT_R_rprime_r = 0.0;
  IT_grad_1_R = 0.0;
  IT_grad_R = 0.0;

  for (int i=0 ; i<3 ; ++i) {
    switch (i)
    {
      case 0: r_plus__i = T.r_nodes (1); r_minus__i = T.r_nodes (0); break;
      case 1: r_plus__i = T.r_nodes (2); r_minus__i = T.r_nodes (1); break;
      case 2: r_plus__i = T.r_nodes (0); r_minus__i = T.r_nodes (2); break;
    }
    // s_plus__i, s_minus__i computation
    const blitz::TinyVector<double, 3> r_plus__i_r(r_plus__i-r), r_minus__i_r(r_minus__i-r);
    s_plus__i = dot(r_plus__i_r, T.s_i_hat (i));
    s_minus__i = dot(r_minus__i_r, T.s_i_hat (i));

    // t_i_0 : distance from r (projected on plane of triangle) to edge
    t_i_0 = dot(r_plus__i_r, T.m_i_hat (i));

    // R_plus__i, R_minus__i, R_i_0 computation
    R_plus__i = sqrt(dot(r_plus__i_r, r_plus__i_r));
    R_minus__i = sqrt(dot(r_minus__i_r, r_minus__i_r));
    R_i_0_square = t_i_0*t_i_0 + h*h;

    // different cases according to the position vector    
    if (abs(t_i_0) > 1.0e-8) {
      beta_i = atan(t_i_0*s_plus__i/(R_i_0_square + abs_h*R_plus__i)) - atan(t_i_0*s_minus__i/(R_i_0_square + abs_h*R_minus__i));
      I_L_minus_1__i = log((R_plus__i+s_plus__i)/(R_minus__i+s_minus__i));
      I_L_plus_1__i = 0.5 * (s_plus__i*R_plus__i - s_minus__i*R_minus__i + R_i_0_square*I_L_minus_1__i);
    }
    else { // if (abs(t_i_0)<1.0e-8)
      t_i_0 = 0.0;
      beta_i = 0.0;
      if (abs_h<1.0e-8) {
        // the following line is derived with the help of HOSPITAL rule
        I_L_minus_1__i = (s_minus__i*s_plus__i <= 0.0) ? 1.0e+90 : s_plus__i/abs(s_plus__i) * log(s_plus__i/s_minus__i);
        I_L_plus_1__i = 0.5 * (s_plus__i*R_plus__i - s_minus__i*R_minus__i);
      }
      else { // if (abs_h>1.0e-8)
        I_L_minus_1__i = log((R_plus__i+s_plus__i)/(R_minus__i+s_minus__i));
        I_L_plus_1__i = 0.5 * (s_plus__i*R_plus__i - s_minus__i*R_minus__i + R_i_0_square*I_L_minus_1__i);
      }
    }
    I_L_plus_3__i = 0.25 * (s_plus__i*R_plus__i*R_plus__i*R_plus__i - s_minus__i*R_minus__i*R_minus__i*R_minus__i + 3.0 * R_i_0_square * I_L_plus_1__i);

    K_1_minus_1__i = t_i_0*I_L_minus_1__i - abs_h * beta_i;
    K_1_plus_1__i = 1.0/3.0 * (h*h * K_1_minus_1__i + t_i_0*I_L_plus_1__i);
    K_2_minus_1__i = T.m_i_hat (i) * I_L_plus_1__i - T.n_hat * (h * K_1_minus_1__i); 
    K_2_plus_1__i = 1.0/3.0 * T.m_i_hat (i) * I_L_plus_3__i - T.n_hat * (h * K_1_plus_1__i);
    K_3_minus_1__i = (-sign_h * beta_i) * T.n_hat - I_L_minus_1__i * T.m_i_hat (i);
    K_3_plus_1__i = (h * K_1_minus_1__i) * T.n_hat - I_L_plus_1__i * T.m_i_hat (i);

    IT_1_R += K_1_minus_1__i;
    IT_R += K_1_plus_1__i;
    IT_1_R_rprime_r += K_2_minus_1__i;
    IT_R_rprime_r += K_2_plus_1__i;
    IT_grad_1_R += K_3_minus_1__i;
    IT_grad_R += K_3_plus_1__i;
  }
}

void ITs_free (std::complex<double>& ITs_G,
               blitz::TinyVector<std::complex<double>, 3>& ITs_G_rprime_r,
               blitz::TinyVector<std::complex<double>, 3>& ITs_grad_G,
               const blitz::TinyVector<double,3>& r,
               const Triangle & Ts,
               const std::complex<double> k,
               const int N_points,
               const int EXTRACT_1_R,
               const int EXTRACT_R)
{
  int j;
  double sum_weigths, norm_factor, R, IT_1_R, IT_R;
  std::complex<double> G_j, I_k_R, exp_minus_I_k_R, k_square;
  blitz::TinyVector<double, 3> rprime, rprime_r, IT_1_R_rprime_r, IT_R_rprime_r, IT_grad_1_R, IT_grad_R;

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points);
  norm_factor = Ts.A/sum_weigths;

  ITs_G = 0.0; // complex<double>
  ITs_G_rprime_r = 0.0; // TinyVector<complex<double>, 3>
  ITs_grad_G = 0.0; // TinyVector<complex<double>, 3>
  if ((EXTRACT_1_R==0) && (EXTRACT_R==0)) { // no singularity extraction
    for (j=0 ; j<N_points ; j++) {
      rprime = Ts.r_nodes (0)*xi[j] + Ts.r_nodes (1)*eta[j] + Ts.r_nodes (2)*(1.0-xi[j]-eta[j]);
      rprime_r = rprime-r;
      R = sqrt(dot(rprime_r, rprime_r));
      I_k_R = I*k*R;
      exp_minus_I_k_R = exp(-I_k_R);
      G_j = exp_minus_I_k_R/R * weigths[j];
      ITs_G += G_j;
      ITs_G_rprime_r += G_j * rprime_r;
      ITs_grad_G += (G_j * (1.0+I_k_R)/(R*R)) * rprime_r;
    }
    ITs_G *= norm_factor;
    ITs_G_rprime_r *= norm_factor;
    ITs_grad_G *= norm_factor;
  }
 
  else if ((EXTRACT_1_R==1) && (EXTRACT_R==0)) { // 1/R singularity extraction
    for (j=0 ; j<N_points ; j++) {
      rprime = Ts.r_nodes (0)*xi[j] + Ts.r_nodes (1)*eta[j] + Ts.r_nodes (2)*(1-xi[j]-eta[j]);
      rprime_r = rprime-r;
      R = sqrt(dot(rprime_r, rprime_r));
      I_k_R = I*k*R;
      exp_minus_I_k_R = exp(-I_k_R);
      G_j = (R>1.0e-10) ? (exp_minus_I_k_R - 1.0)/R * weigths[j] : -I * k * weigths[j];
      ITs_G += G_j;
      ITs_G_rprime_r += G_j * rprime_r;
      if (R>1.0e-10) ITs_grad_G -= ( (-exp_minus_I_k_R*(1.0+I_k_R) + 1.0)/(R*R*R) * weigths[j] ) * rprime_r;
    }
    IT_singularities (IT_1_R, IT_R, IT_1_R_rprime_r, IT_R_rprime_r, IT_grad_1_R, IT_grad_R, r, Ts);
    ITs_G = ITs_G * norm_factor + IT_1_R;
    ITs_G_rprime_r = ITs_G_rprime_r * norm_factor + IT_1_R_rprime_r;
    ITs_grad_G = ITs_grad_G * norm_factor + IT_grad_1_R;
  }

  else if ((EXTRACT_1_R==1) && (EXTRACT_R==1)) { // 1/R and R singularity extraction
    k_square = k*k;
    for (j=0 ; j<N_points ; j++) {
      rprime = Ts.r_nodes (0)*xi[j] + Ts.r_nodes (1)*eta[j] + Ts.r_nodes (2)*(1-xi[j]-eta[j]);
      rprime_r = rprime-r;
      R = sqrt(dot(rprime_r, rprime_r));
      I_k_R = I*k*R;
      exp_minus_I_k_R = exp(-I_k_R);
      G_j = (R>1.0e-10) ? ( (exp_minus_I_k_R - 1.0)/R + k_square/2.0 * R ) * weigths[j] : -I * k * weigths[j];
      ITs_G += G_j;
      ITs_G_rprime_r += G_j * rprime_r;
      if (R>1.0e-10) ITs_grad_G -= ( (-exp_minus_I_k_R*(1.0+I_k_R) + 1.0 + k_square/2.0 * R*R)/(R*R*R) * weigths[j] ) * rprime_r;
    }
    IT_singularities (IT_1_R, IT_R, IT_1_R_rprime_r, IT_R_rprime_r, IT_grad_1_R, IT_grad_R, r, Ts);
    ITs_G = ITs_G * norm_factor + IT_1_R - k_square/2.0 * IT_R;
    ITs_G_rprime_r = ITs_G_rprime_r * norm_factor + IT_1_R_rprime_r - k_square/2.0 * IT_R_rprime_r;
    ITs_grad_G = ITs_grad_G * norm_factor + IT_grad_1_R - k_square/2.0 * IT_grad_R;
  }
}

void ITo_ITs_free (std::complex<double>& ITo_ITs_G,
                   blitz::TinyVector<std::complex<double>, 3>& ITo_r_ITs_G,
                   blitz::TinyVector<std::complex<double>, 3>& ITo_ITs_G_rprime,
                   std::complex<double>& ITo_r_dot_ITs_G_rprime,
                   blitz::TinyVector<std::complex<double>, 3>& ITo_n_hat_X_r_ITs_G,
                   std::complex<double>& ITo_n_hat_X_r_dot_ITs_G_rprime,
                   blitz::TinyVector<std::complex<double>, 3>& ITo_ITs_grad_G,
                   blitz::TinyVector<std::complex<double>, 3>& ITo_r_X_ITs_grad_G,
                   std::complex<double> & ITo_n_hat_X_r_dot_r_X_ITs_grad_G,
                   blitz::TinyVector<std::complex<double>, 3>& ITo_n_hat_X_r_X_ITs_grad_G,
                   const Triangle & To,
                   const Triangle & Ts,
                   const std::complex<double> k,
                   const int N_points_o,
                   const int N_points_s,
                   const int EXTRACT_1_R,
                   const int EXTRACT_R)
{
  int j;
  double sum_weigths, norm_factor;
  std::complex<double> ITs_G_j;
  blitz::TinyVector<double, 3> r, n_hat_X_r;
  blitz::TinyVector<std::complex<double>, 3> r_ITs_G_j, ITs_G_rprime_r_j, ITs_G_rprime_j, ITs_grad_G_j, r_X_ITs_grad_G_j, n_hat_X_r_X_ITs_grad_G_j;

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points_o);

  ITo_ITs_G = 0.0; // complex<double>
  ITo_r_ITs_G = 0.0; // TinyVector<complex<double>, 3>
  ITo_ITs_G_rprime = 0.0; // TinyVector<complex<double>, 3>
  ITo_r_dot_ITs_G_rprime = 0.0; // complex<double>
  ITo_n_hat_X_r_ITs_G = 0.0; // TinyVector<complex<double>, 3>
  ITo_n_hat_X_r_dot_ITs_G_rprime = 0.0; // complex<double>
  ITo_ITs_grad_G = 0.0; // TinyVector<complex<double>, 3>
  ITo_r_X_ITs_grad_G = 0.0; // TinyVector<complex<double>, 3>
  ITo_n_hat_X_r_dot_r_X_ITs_grad_G = 0.0; // complex<double>
  ITo_n_hat_X_r_X_ITs_grad_G = 0.0; // TinyVector<complex<double>, 3>
  for (j=0 ; j<N_points_o ; j++) {
    r = To.r_nodes (0)*xi[j] + To.r_nodes (1)*eta[j] + To.r_nodes (2)*(1-xi[j]-eta[j]);
    n_hat_X_r = cross(To.n_hat, r);
    ITs_free (ITs_G_j, ITs_G_rprime_r_j, ITs_grad_G_j, r, Ts, k, N_points_s, EXTRACT_1_R, EXTRACT_R);

    ITs_G_j *= weigths[j];
    r_ITs_G_j = ITs_G_j * r;
    ITs_G_rprime_j = ITs_G_rprime_r_j * weigths[j] + r_ITs_G_j;
    ITs_grad_G_j *= weigths[j];
    r_X_ITs_grad_G_j = r(1) * ITs_grad_G_j(2) - r(2) * ITs_grad_G_j(1),
                       r(2) * ITs_grad_G_j(0) - r(0) * ITs_grad_G_j(2),
                       r(0) * ITs_grad_G_j(1) - r(1) * ITs_grad_G_j(0);
    n_hat_X_r_X_ITs_grad_G_j = n_hat_X_r(1) * ITs_grad_G_j(2) - n_hat_X_r(2) * ITs_grad_G_j(1),
                               n_hat_X_r(2) * ITs_grad_G_j(0) - n_hat_X_r(0) * ITs_grad_G_j(2),
                               n_hat_X_r(0) * ITs_grad_G_j(1) - n_hat_X_r(1) * ITs_grad_G_j(0);

    ITo_ITs_G += ITs_G_j;
    ITo_r_ITs_G += r_ITs_G_j;
    ITo_ITs_G_rprime += ITs_G_rprime_j;
    ITo_r_dot_ITs_G_rprime += dot(r, ITs_G_rprime_j);
    ITo_n_hat_X_r_ITs_G += ITs_G_j * n_hat_X_r;
    ITo_n_hat_X_r_dot_ITs_G_rprime += dot(n_hat_X_r, ITs_G_rprime_j);
    ITo_ITs_grad_G += ITs_grad_G_j;
    ITo_r_X_ITs_grad_G +=  r_X_ITs_grad_G_j;
    ITo_n_hat_X_r_dot_r_X_ITs_grad_G += dot (n_hat_X_r, r_X_ITs_grad_G_j);
    ITo_n_hat_X_r_X_ITs_grad_G += n_hat_X_r_X_ITs_grad_G_j;
  }
  norm_factor = To.A/(4.0*M_PI*sum_weigths);
  ITo_ITs_G *= norm_factor;
  ITo_r_ITs_G *= norm_factor;
  ITo_ITs_G_rprime *= norm_factor;
  ITo_r_dot_ITs_G_rprime *= norm_factor;
  ITo_n_hat_X_r_ITs_G *= norm_factor;
  ITo_n_hat_X_r_dot_ITs_G_rprime *= norm_factor;
  ITo_ITs_grad_G *= norm_factor;
  ITo_r_X_ITs_grad_G *= norm_factor;
  ITo_n_hat_X_r_dot_r_X_ITs_grad_G *= norm_factor;
  ITo_n_hat_X_r_X_ITs_grad_G *= norm_factor;
}

void IDTo_ITs_free (std::complex<double> & IDTo_l_hat_dot_r_ITs_G,
                    blitz::TinyVector<complex<double>, 3>& IDTo_l_hat_ITs_G,
                    const Triangle & To,
                    const Triangle & Ts,
                    const std::complex<double> k,
                    const int N_points_o,
                    const int N_points_s,
                    const int EXTRACT_1_R,
                    const int EXTRACT_R)
{
  blitz::Array<double, 1> XGL, WGL;
  Gauss_Legendre (XGL, WGL, N_points_o);

  double norm_factor, norm_r_hlgth;
  blitz::TinyVector<double, 3> r, r_center, n_hat_X_r, m_hat, l_hat, r_plus__i, r_minus__i, r_hlgth;

  std::complex<double> ITs_G_j, I_k;
  blitz::TinyVector<std::complex<double>, 3> ITs_G_rprime_r_j, ITs_grad_G_j, I_r_k;

  IDTo_l_hat_dot_r_ITs_G = 0.0; // complex<double>
  IDTo_l_hat_ITs_G = 0.0; // TinyVector<complex<double>, 3>
  for (int i=0 ; i<3 ; i++) { // we turn on the contour of T_obs
    switch (i)
    {
      case 0: r_plus__i = To.r_nodes (1); r_minus__i = To.r_nodes (0); break;
      case 1: r_plus__i = To.r_nodes (2); r_minus__i = To.r_nodes (1); break;
      case 2: r_plus__i = To.r_nodes (0); r_minus__i = To.r_nodes (2); break;
    }
    l_hat = To.s_i_hat (i);

    // This is a Gauss-Kronrod rule, applied to each edge
    r_hlgth = 0.5*(r_plus__i-r_minus__i);
    norm_r_hlgth = sqrt(dot(r_hlgth, r_hlgth));
    r_center = 0.5*(r_plus__i+r_minus__i);
    I_k = 0.0;
    I_r_k = 0.0;
    for (int j=0 ; j<N_points_o ; j++) {
      r = r_center + r_hlgth * XGL(j);
      ITs_free (ITs_G_j, ITs_G_rprime_r_j, ITs_grad_G_j, r, Ts, k, N_points_s, EXTRACT_1_R, EXTRACT_R);
      I_k += ITs_G_j * WGL(j);
      I_r_k += (ITs_G_j * WGL(j)) * r;
    } // end of the Gauss-Kronrod rule applied to edge (i)
    IDTo_l_hat_ITs_G += I_k * norm_r_hlgth * l_hat;
    IDTo_l_hat_dot_r_ITs_G += dot (l_hat, I_r_k) * norm_r_hlgth;
  } // contour integration finished
  norm_factor = 1.0/(4.0*M_PI);
  IDTo_l_hat_dot_r_ITs_G *= norm_factor;
  IDTo_l_hat_ITs_G *= norm_factor;
}

/* special functions for the excitation vectors calculations */
void V_EH_ITo_free (std::complex<double>& ITo_G,
                    blitz::TinyVector<std::complex<double>, 3>& ITo_G_rprime_r,
                    blitz::TinyVector<std::complex<double>, 3>& ITo_grad_G,
                    blitz::TinyVector<std::complex<double>, 3>& ITo_n_hat_X_r_X_grad_G,
                    const blitz::TinyVector<double,3>& r,
                    const Triangle & To,
                    const std::complex<double> k,
                    const int N_points,
                    const int EXTRACT_1_R,
                    const int EXTRACT_R)
{
  int j;
  double sum_weigths, norm_factor, R, IT_1_R, IT_R;
  std::complex<double> G_j, I_k_R, exp_minus_I_k_R, k_square;
  blitz::TinyVector<double, 3> rprime, rprime_r, IT_1_R_rprime_r, IT_R_rprime_r, IT_grad_1_R, IT_grad_R, n_hat_X_rprime;
  blitz::TinyVector<std::complex<double>, 3> ITo_grad_G_j;

  const double *xi, *eta, *weigths;
  IT_points (xi, eta, weigths, sum_weigths, N_points);
  norm_factor = To.A/sum_weigths;

  ITo_G = 0.0; // complex<double>
  ITo_G_rprime_r = 0.0; // TinyVector<complex<double>, 3>
  ITo_grad_G = 0.0; // TinyVector<complex<double>, 3>
  ITo_n_hat_X_r_X_grad_G = 0.0; // TinyVector<complex<double>, 3>
  if ((EXTRACT_1_R==0) && (EXTRACT_R==0)) { // no singularity extraction
    for (j=0 ; j<N_points ; j++) {
      rprime = To.r_nodes (0)*xi[j] + To.r_nodes (1)*eta[j] + To.r_nodes (2)*(1.0-xi[j]-eta[j]);
      rprime_r = rprime-r;
      R = sqrt(dot(rprime_r, rprime_r));
      I_k_R = I*k*R;
      exp_minus_I_k_R = exp(-I_k_R);
      G_j = exp_minus_I_k_R/R * weigths[j];
      ITo_G += G_j;
      ITo_G_rprime_r += G_j * rprime_r;
      ITo_grad_G_j = (G_j * (1.0+I_k_R)/(R*R)) * rprime_r;
      ITo_grad_G += ITo_grad_G_j;
      n_hat_X_rprime = To.n_hat(1)*rprime(2)-To.n_hat(2)*rprime(1),
                       To.n_hat(2)*rprime(0)-To.n_hat(0)*rprime(2),
                       To.n_hat(0)*rprime(1)-To.n_hat(1)*rprime(0);
      ITo_n_hat_X_r_X_grad_G += n_hat_X_rprime(1) * ITo_grad_G_j(2) - n_hat_X_rprime(2) * ITo_grad_G_j(1),
                                n_hat_X_rprime(2) * ITo_grad_G_j(0) - n_hat_X_rprime(0) * ITo_grad_G_j(2),
                                n_hat_X_rprime(0) * ITo_grad_G_j(1) - n_hat_X_rprime(1) * ITo_grad_G_j(0);
    }
    ITo_G *= norm_factor;
    ITo_G_rprime_r *= norm_factor;
    ITo_grad_G *= norm_factor;
    ITo_n_hat_X_r_X_grad_G *= norm_factor;
  }
 
  else if ((EXTRACT_1_R==1) && (EXTRACT_R==0)) { // 1/R singularity extraction
    for (j=0 ; j<N_points ; j++) {
      rprime = To.r_nodes (0)*xi[j] + To.r_nodes (1)*eta[j] + To.r_nodes (2)*(1-xi[j]-eta[j]);
      rprime_r = rprime-r;
      R = sqrt(dot(rprime_r, rprime_r));
      I_k_R = I*k*R;
      exp_minus_I_k_R = exp(-I_k_R);
      G_j = (R>1.0e-10) ? (exp_minus_I_k_R - 1.0)/R * weigths[j] : -I * k * weigths[j];
      ITo_G += G_j;
      ITo_G_rprime_r += G_j * rprime_r;
      if (R>1.0e-10) ITo_grad_G_j = ( (-exp_minus_I_k_R*(1.0+I_k_R) + 1.0)/(R*R*R) * weigths[j] ) * rprime_r;
      ITo_grad_G -= ITo_grad_G_j;
      n_hat_X_rprime = To.n_hat(1)*rprime(2)-To.n_hat(2)*rprime(1),
                       To.n_hat(2)*rprime(0)-To.n_hat(0)*rprime(2),
                       To.n_hat(0)*rprime(1)-To.n_hat(1)*rprime(0);
      ITo_n_hat_X_r_X_grad_G += n_hat_X_rprime(1) * ITo_grad_G_j(2) - n_hat_X_rprime(2) * ITo_grad_G_j(1),
                                n_hat_X_rprime(2) * ITo_grad_G_j(0) - n_hat_X_rprime(0) * ITo_grad_G_j(2),
                                n_hat_X_rprime(0) * ITo_grad_G_j(1) - n_hat_X_rprime(1) * ITo_grad_G_j(0);
    }
    IT_singularities (IT_1_R, IT_R, IT_1_R_rprime_r, IT_R_rprime_r, IT_grad_1_R, IT_grad_R, r, To);
    ITo_G = ITo_G * norm_factor + IT_1_R;
    ITo_G_rprime_r = ITo_G_rprime_r * norm_factor + IT_1_R_rprime_r;
    ITo_grad_G = ITo_grad_G * norm_factor + IT_grad_1_R;
    ITo_n_hat_X_r_X_grad_G = ITo_n_hat_X_r_X_grad_G * norm_factor;
  }

  else if ((EXTRACT_1_R==1) && (EXTRACT_R==1)) { // 1/R and R singularity extraction
    k_square = k*k;
    for (j=0 ; j<N_points ; j++) {
      rprime = To.r_nodes (0)*xi[j] + To.r_nodes (1)*eta[j] + To.r_nodes (2)*(1-xi[j]-eta[j]);
      rprime_r = rprime-r;
      R = sqrt(dot(rprime_r, rprime_r));
      I_k_R = I*k*R;
      exp_minus_I_k_R = exp(-I_k_R);
      G_j = (R>1.0e-10) ? ( (exp_minus_I_k_R - 1.0)/R + k_square/2.0 * R ) * weigths[j] : -I * k * weigths[j];
      ITo_G += G_j;
      ITo_G_rprime_r += G_j * rprime_r;
      if (R>1.0e-10) ITo_grad_G_j = ( (-exp_minus_I_k_R*(1.0+I_k_R) + 1.0 + k_square/2.0 * R*R)/(R*R*R) * weigths[j] ) * rprime_r;
      ITo_grad_G -= ITo_grad_G_j;
      n_hat_X_rprime = To.n_hat(1)*rprime(2)-To.n_hat(2)*rprime(1),
                       To.n_hat(2)*rprime(0)-To.n_hat(0)*rprime(2),
                       To.n_hat(0)*rprime(1)-To.n_hat(1)*rprime(0);
      ITo_n_hat_X_r_X_grad_G += n_hat_X_rprime(1) * ITo_grad_G_j(2) - n_hat_X_rprime(2) * ITo_grad_G_j(1),
                                n_hat_X_rprime(2) * ITo_grad_G_j(0) - n_hat_X_rprime(0) * ITo_grad_G_j(2),
                                n_hat_X_rprime(0) * ITo_grad_G_j(1) - n_hat_X_rprime(1) * ITo_grad_G_j(0);
    }
    IT_singularities (IT_1_R, IT_R, IT_1_R_rprime_r, IT_R_rprime_r, IT_grad_1_R, IT_grad_R, r, To);
    ITo_G = ITo_G * norm_factor + IT_1_R - k_square/2.0 * IT_R;
    ITo_G_rprime_r = ITo_G_rprime_r * norm_factor + IT_1_R_rprime_r - k_square/2.0 * IT_R_rprime_r;
    ITo_grad_G = ITo_grad_G * norm_factor + IT_grad_1_R - k_square/2.0 * IT_grad_R;
    //ITo_n_hat_X_r_X_grad_G = ITo_n_hat_X_r_X_grad_G * norm_factor;
  }
}


