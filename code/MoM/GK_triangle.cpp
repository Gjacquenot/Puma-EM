using namespace std;

#include "GK_triangle.h"

const int MAX_N_POINTS = 9;
const double xi_1[1] = {0.333333333333333};
const double eta_1[1] = {0.333333333333333};
const double weigths_1[1] = {1.0};
const double sum_weigths_1 = 1.0;

const double xi_3[3] = { 0.0, 0.5, 0.5 };
const double eta_3[3] = { 0.5, 0.0, 0.5 };
const double weigths_3[3] = { 0.333333333333333, 0.333333333333333, 0.333333333333333 };
const double sum_weigths_3 = 1.0;

const double xi_6[6] = { 0.816847572980459, 0.091576213509771, 0.091576213509771, 0.108103018168070, 0.445948490915965, 0.445948490915965 };
const double eta_6[6] = { 0.091576213509771, 0.816847572980459, 0.091576213509771, 0.445948490915965, 0.108103018168070, 0.445948490915965 };
const double weigths_6[6] = { 0.109951743655322, 0.109951743655322, 0.109951743655322, 0.223381589678011, 0.223381589678011, 0.223381589678011 };
const double sum_weigths_6 = 1.0;

const double xi_9[9] = { 0.124949503233232, 0.437525248383384, 0.437525248383384, 0.797112651860071, 0.797112651860071, 0.165409927389841, 0.165409927389841, 0.037477420750088, 0.037477420750088 };
const double eta_9[9] = { 0.437525248383384, 0.124949503233232, 0.437525248383384, 0.165409927389841, 0.037477420750088, 0.797112651860071, 0.037477420750088, 0.797112651860071, 0.165409927389841 };
const double weigths_9[9] = { 0.205950504760887, 0.205950504760887, 0.205950504760887, 0.063691414286223, 0.063691414286223, 0.063691414286223, 0.063691414286223, 0.063691414286223, 0.063691414286223 };
const double sum_weigths_9 = 1.0;

const double w1=0.050844906370207, x1 = 0.873821971016996, y1 = 0.063089014491502, w2 = 0.116786275726379, x2 = 0.501426509658179, y2 = 0.249286745170910, w3 = 0.082851075618374, x3 = 0.636502499121399, y3 = 0.310352451033785, z3 = 0.053145049844816;
const double xi_12[12] = { x1, y1, y1, x2, y2, y2, x3, x3, y3, y3, z3, z3 };
const double eta_12[12] = { y1, x1, y1, y2, x2, y2, y3, z3, x3, z3, x3, y3 };
const double weigths_12[12] = { w1, w1, w1, w2, w2, w2, w3, w3, w3, w3, w3, w3 };
const double sum_weigths_12 = 1.0;

const double a = 0.479308067841923, b = 0.260345966079038, c = 0.869739794195568, d = 0.065130102902216, e = 0.638444188569809, f = 0.312865496004875;
const double g = 0.048690315425316, h = 0.333333333333333, t = 0.175615257433204, u = 0.053347235608839, v = 0.077113760890257, w = -0.14957004446767;
const double xi_13[13] = {a, b, b, c, d, d, e, e, f, f, g, g, h};
const double eta_13[13] = {b, a, b, d, c, d, f, g, e, g, e, f, h};
const double weigths_13[13] = {t, t, t, u, u, u, v, v, v, v, v, v, w};
const double sum_weigths_13 = 1.0; 

void IT_points (const double * & xi, const double * & eta, const double * & weigths, double & sum_weigths, const int N) {

  switch (N)
  {
    case 1: xi = xi_1; eta = eta_1; weigths = weigths_1; sum_weigths = sum_weigths_1; break;
    case 3: xi = xi_3; eta = eta_3; weigths = weigths_3; sum_weigths = sum_weigths_3; break;
    case 6: xi = xi_6; eta = eta_6; weigths = weigths_6; sum_weigths = sum_weigths_6; break;
    case 9: xi = xi_9; eta = eta_9; weigths = weigths_9; sum_weigths = sum_weigths_9; break;
    case 12: xi = xi_12; eta = eta_12; weigths = weigths_12; sum_weigths = sum_weigths_12; break;
    case 13: xi = xi_13; eta = eta_13; weigths = weigths_13; sum_weigths = sum_weigths_13; break;
  }
}
