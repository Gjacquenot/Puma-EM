#ifndef TRIANGLE_INT_H
#define TRIANGLE_INT_H
#include <vector>

using namespace std;

#include "dictionary.h"

template <class T>
inline void cross3D(T result[], const T a[], const T b[])
{
  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
}

template <class T>
inline T dot3D(const T a[], const T b[])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

class Triangle {
  public:
    int number;
    std::vector<int> RWGIndexes;
    std::vector<int> indexesInRWGs;
    std::vector<double> signInRWG;
    // triangle's nodes
    double r_nodes[3][3];
    // triangle's normals to edges
    double m_i_hat[3][3];
    // triangle's edges unit vectors
    double s_i_hat[3][3];
    // triangle's normal unit vector
    double n_hat[3];
    // triangle's gravity centre
    double r_grav[3];
    double A;
    double R_max;

    // constructors    
    Triangle(void){};
    Triangle(const double r0[],
             const double r1[],
             const double r2[],
             const int tr_number);
    //! allows a triangle to be overwritten by the provided triangle
    /*!
      \param triangleToCopy the triangle to copy
    */
    void copyTriangle (const Triangle& triangleToCopy);
    Triangle(const Triangle&); // copy constructor
    Triangle& operator=(const Triangle&); // copy assignment operator
    //! the destructor
    ~Triangle();
};

class RWG {
  public:
    int number;
    int triangleNumbers[2];
    //! the signs of the triangles for the basis function
    int triangleSigns[2];
    double length;
    //! the coordinates of the opposite vector to the RWG edge in the given triangle
    double vertexesCoord_0[3];
    double vertexesCoord_1[3];
    double vertexesCoord_2[3];
    double vertexesCoord_3[3];
    // constructors
    RWG(void){};
    RWG(const int RWG_number,
        const int triangle_numbers[],
        const int triangle_signs[],
        const double r0[],
        const double r1[],
        const double r2[],
        const double r3[]);
    //! allows a RWG_function to be overwritten by the provided RWG_function
    /*!
      \param RWG_functionToCopy the RWG_function to copy
    */
    void copyRWG (const RWG & RWGToCopy);
    RWG(const RWG & RWGToCopy); // copy constructor
    RWG& operator=(const RWG & RWGToCopy); // copy assignment operator
    //! the destructor
    ~RWG(){};
};

void constructVectorTriangles(std::vector<Triangle>& triangles,
                              const std::vector<RWG>& vectorRWGs,
                              const std::vector<Dictionary<int, int> >& TriangleToRWG);

void IT_fm_fn (double & IT_r_square, 
               double IT_r[], // dim 3
               const Triangle & T);

void ITs_free (std::complex<double>& ITs_G,
               std::complex<double> ITs_G_rprime_r[], // dim 3
               std::complex<double> ITs_grad_G[], // dim 3
               const double r[], // dim 3
               const Triangle & Ts,
               const std::complex<double> & k,
               const int N_points,
               const double * xi,
               const double * eta,
               const double * weights,
               const double sum_weights,
               const int EXTRACT_1_R,
               const int EXTRACT_R);

void ITo_ITs_free (std::complex<double>& ITo_ITs_G, 
                   std::complex<double> ITo_r_ITs_G[], // dim 3 
                   std::complex<double> ITo_ITs_G_rprime[], // dim 3
                   std::complex<double>& ITo_r_dot_ITs_G_rprime, 
                   std::complex<double> ITo_n_hat_X_r_ITs_G[], // dim 3
                   std::complex<double>& ITo_n_hat_X_r_dot_ITs_G_rprime, 
                   std::complex<double> ITo_ITs_grad_G[], // dim 3
                   std::complex<double> ITo_r_X_ITs_grad_G[], // dim 3
                   std::complex<double> & ITo_n_hat_X_r_dot_r_X_ITs_grad_G, 
                   std::complex<double> ITo_n_hat_X_r_X_ITs_grad_G[], // dim 3
                   const Triangle & To, 
                   const Triangle & Ts, 
                   const std::complex<double> k, 
                   const int N_points_o, 
                   const int N_points_s, 
                   const int EXTRACT_1_R, 
                   const int EXTRACT_R);

void IDTo_ITs_free (std::complex<double> & IDTo_l_hat_dot_r_ITs_G, 
                    std::complex<double> IDTo_l_hat_ITs_G[], // dim 3
                    const Triangle & To, 
                    const Triangle & Ts, 
                    const std::complex<double> k, 
                    const int N_points_o, 
                    const int N_points_s, 
                    const int EXTRACT_1_R, 
                    const int EXTRACT_R);

/* special functions for the excitation vectors calculations */
void V_EH_ITo_free (std::complex<double>& ITo_G,
                    std::complex<double> ITo_G_rprime_r[], // dim 3
                    std::complex<double> ITo_grad_G[], // dim 3
                    std::complex<double> ITo_n_hat_X_r_X_grad_G[], // dim 3
                    double r[], // dim 3
                    const Triangle & To,
                    const std::complex<double> k,
                    const int N_points,
                    const int EXTRACT_1_R,
                    const int EXTRACT_R);

#endif
