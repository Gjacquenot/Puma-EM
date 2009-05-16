/***************************************************************************
 * GK.cpp  Gauss-Kronrod integration routine specialized for spectral DGFs
 *         it is used in the Sommerfeld integrals
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
#include <float.h>
#include <complex>
#include <blitz/array.h>

using namespace blitz;

const complex<double> I (0.0, 1.0);
const double XGK15[8] = {
		0.99145537112081263921,
		0.94910791234275852453,
		0.86486442335976907279,
		0.74153118559939443986,
		0.58608723546769113029,
		0.40584515137739716691,
		0.20778495500789846760,
		0.00000000000000000000};
const double WGK15[8] = {
		0.02293532201052922496,
		0.06309209262997855329,
		0.10479001032225018384,
		0.14065325971552591875,
		0.16900472663926790283,
		0.19035057806478540991,
		0.20443294007529889241,
		0.20948214108472782801};
const double WG7[4] = {
		0.12948496616886969327,
		0.27970539148927666790,
		0.38183005050511894495,
		0.41795918367346938776};
const double XGK31[16] = {
		0.99800229869339706029,
		0.98799251802048542849,
		0.96773907567913913426,
		0.93727339240070590431,
		0.89726453234408190088,
		0.84820658341042721620,
		0.79041850144246593297,
		0.72441773136017004742,
		0.65099674129741697053,
		0.57097217260853884754,
		0.48508186364023968069,
		0.39415134707756336990,
		0.29918000715316881217,
		0.20119409399743452230,
		0.10114206691871749903,
		0.00000000000000000000};
const double WGK31[16] = {
		0.00537747987292334899,
		0.01500794732931612254,
		0.02546084732671532019,
		0.03534636079137584622,
		0.04458975132476487661,
		0.05348152469092808727,
		0.06200956780067064029,
		0.06985412131872825871,
		0.07684968075772037889,
		0.08308050282313302104,
		0.08856444305621177065,
		0.09312659817082532123,
		0.09664272698362367851,
		0.09917359872179195933,
		0.10076984552387559504,
		0.10133000701479154902};
const double WG15[8] = {
		0.03075324199611726835,
		0.07036604748810812471,
		0.10715922046717193501,
		0.13957067792615431445,
		0.16626920581699393355,
		0.18616100001556221103,
		0.19843148532711157646,
		0.20257824192556127288};

#include "layers_constants.h"

complex<double> GK15 (complex<double> (*f)(const complex<double>, const double, const double, const double, const double, const int, const int, const layers_constants &), const double v, const complex<double> a, const complex<double> b, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> hlgth;
  complex<double> absc, centr;
  complex<double> fc, fsum, fval1, fval2;
  complex<double> resg, resk, reskh, result;
  int j, jtw, jtwm1;
  
  centr = 0.5 * (a+b);
  hlgth = 0.5 * (b-a);
 
  fc = f (centr, v, rho, z, z_prime, m, n, LC);
  resg = fc * WG7[3];
  resk = fc * WGK15[7];
  for (j=0 ; j<3 ; j++) {
    jtw = 2*j + 1;
    absc = hlgth * XGK15[jtw];
    fval1 = f (centr-absc, v, rho, z, z_prime, m, n, LC);
    fval2 = f (centr+absc, v, rho, z, z_prime, m, n, LC);
    fsum = fval1 + fval2;
    resg = resg + WG7[j] * fsum;
    resk = resk + WGK15[jtw] * fsum;
  }
  for (j=0 ; j<4 ; j++) { 
    jtwm1 = j*2;
    absc = hlgth * XGK15[jtwm1];
    fval1 = f (centr-absc, v, rho, z, z_prime, m, n, LC);
    fval2 = f (centr+absc, v, rho, z, z_prime, m, n, LC);
    fsum = fval1 + fval2;
    resk = resk + WGK15[jtwm1] * fsum;
  }
  reskh = resk * 0.5;
  result = resk * abs(hlgth);
  return result;
}

complex<double> GK31 (complex<double> (*f)(const complex<double>, const double, const double, const double, const double, const int, const int, const layers_constants &), const double v, const complex<double> a, const complex<double> b, const double rho, const double z, const double z_prime, const int m, const int n, const layers_constants & LC) {

  complex<double> hlgth;
  complex<double> absc, centr;
  complex<double> fc, fsum, fval1, fval2;
  complex<double> resg, resk, reskh, result;
  int j, jtw, jtwm1;
  
  centr = 0.5 * (a+b);
  hlgth = 0.5 * (b-a);

  fc = f (centr, v, rho, z, z_prime, m, n, LC);
  resg = fc * WG15[7];
  resk = fc * WGK31[15];
  for (j=0 ; j<7 ; j++) {
    jtw = 2*j + 1;
    absc = hlgth * XGK31[jtw];
    fval1 = f (centr-absc, v, rho, z, z_prime, m, n, LC);
    fval2 = f (centr+absc, v, rho, z, z_prime, m, n, LC);
    fsum = fval1 + fval2;
    resg = resg + WG15[j] * fsum;
    resk = resk + WGK31[jtw] * fsum;
  }
  for (j=0 ; j<8 ; j++) { 
    jtwm1 = j*2;
    absc = hlgth * XGK31[jtwm1];
    fval1 = f (centr-absc, v, rho, z, z_prime, m, n, LC);
    fval2 = f (centr+absc, v, rho, z, z_prime, m, n, LC);
    fsum = fval1 + fval2;
    resk = resk + WGK31[jtwm1] * fsum;
  }
  reskh = resk * 0.5;
  result = resk * abs(hlgth);
  return result;
}
