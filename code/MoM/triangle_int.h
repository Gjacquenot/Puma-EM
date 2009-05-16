#ifndef TRIANGLE_INT_H
#define TRIANGLE_INT_H
#include <vector>
using namespace blitz;

#include "mesh.h"

class Triangle {
  public:
    int number;
    std::vector<int> RWGIndexes;
    std::vector<int> indexesInRWGs;
    std::vector<double> signInRWG;
    blitz::Array<blitz::TinyVector<double, 3>, 1> r_nodes;
    blitz::Array<blitz::TinyVector<double, 3>, 1> m_i_hat;
    blitz::Array<blitz::TinyVector<double, 3>, 1> s_i_hat;
    blitz::TinyVector<double, 3> n_hat;
    blitz::TinyVector<double, 3> r_grav;
    double A;
    double R_max;

    // constructors    
    Triangle(void);
    Triangle(const blitz::TinyVector<double, 3>& r0,
             const blitz::TinyVector<double, 3>& r1,
             const blitz::TinyVector<double, 3>& r2,
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
    blitz::TinyVector<int, 2> triangleNumbers;
    //! the signs of the triangles for the basis function
    blitz::TinyVector<int, 2> triangleSigns;
    double length;
    //! the coordinates of the opposite vector to the RWG edge in the given triangle
    blitz::Array<blitz::TinyVector<double, 3>, 1> vertexesCoord;
    // constructors
    RWG(void){};
    RWG(const int RWG_number,
        const blitz::TinyVector<int, 2>& triangle_numbers,
        const blitz::TinyVector<int, 2>& triangle_signs,
        const blitz::Array<double, 1>& r0,
        const blitz::Array<double, 1>& r1,
        const blitz::Array<double, 1>& r2,
        const blitz::Array<double, 1>& r3);
    //! allows a RWG_function to be overwritten by the provided RWG_function
    /*!
      \param RWG_functionToCopy the RWG_function to copy
    */
    void copyRWG (const RWG & RWGToCopy);
    RWG(const RWG & RWGToCopy); // copy constructor
    RWG& operator=(const RWG & RWGToCopy); // copy assignment operator
    //! the destructor
    ~RWG(){vertexesCoord.free();};
};

void constructVectorTriangles(std::vector<Triangle>& triangles,
                              const std::vector<RWG>& vectorRWGs,
                              const std::vector<Dictionary<int, int> >& TriangleToRWG);

inline void r_node(TinyVector<double,3>& r, const Array<double,2>& vertexes_coord, const int node_index) {
  for (int i=0 ; i<3 ; i++) r(i) = vertexes_coord(node_index, i);
}

void IT_fm_fn (double & IT_r_square, 
               blitz::TinyVector<double, 3>& IT_r, 
               const Triangle & T);

void ITs_free (std::complex<double>& ITs_G,
               blitz::TinyVector<std::complex<double>, 3>& ITs_G_rprime_r,
               blitz::TinyVector<std::complex<double>, 3>& ITs_grad_G,
               const blitz::TinyVector<double,3>& r,
               const Triangle & Ts,
               const std::complex<double> k,
               const int N_points,
               const int EXTRACT_1_R,
               const int EXTRACT_R);

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
                   const int EXTRACT_R);

void IDTo_ITs_free (std::complex<double> & IDTo_l_hat_dot_r_ITs_G, 
                    blitz::TinyVector<std::complex<double>, 3>& IDTo_l_hat_ITs_G, 
                    const Triangle & To, 
                    const Triangle & Ts, 
                    const std::complex<double> k, 
                    const int N_points_o, 
                    const int N_points_s, 
                    const int EXTRACT_1_R, 
                    const int EXTRACT_R);

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
                    const int EXTRACT_R);

#endif
