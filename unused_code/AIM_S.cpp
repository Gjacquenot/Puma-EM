/***************************************************************************
 * AIM_S.cpp  computes the inverse of the Vandermonde matrix
 *            Based upon Bleszynski paper, Radio Science September-October 1996
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
#include <complex>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

void AIM_S (Array<double, 2>& S_x, const Array<double, 1>& x) {

  /* this function computes the inverse of the Vandermonde matrix:
   *
   * [x.^0 ; x.^1 ; ... ; x.^M]
   *
   *  where x is a horizontal vector of positions along a line of size (M+1)
   */

  int i, j, l, M = x.size() - 1, LINES, COLS, line, col, ind;
  double Num, Den, Sign;
  Array<int, 1> k1_k2__kl, possible_indexes_for_x (M);
  Array<int, 2> K;
  Array<double, 1> Num_tmp (M), Den_tmp (M);
  
  //S_x.resize (M+1, M+1);  

  for (j=0 ; j<M+1 ; j++) { // line index
    Sign = -1.0;
    for (l=0 ; l<M+1 ; l++) { // M - column index

      // numerator
      Sign *= -1.0; // Sign = (-1.0)^l
      if (l==0) Num = 1.0;

      else if ( (l==1) || (l==M) ) {
        for (i=0 ; i<M+1 ; i++) { // we "jump" above x (j)
          if (i<j) Num_tmp (i) = x (i);
          else if (i>j) Num_tmp (i-1) = x (i);
	}
        if (l==1) Num = sum (Num_tmp);
        else Num = product (Num_tmp);
      } // (l==1) or (l==M)

      else { // (l!=0, 1, M): the tricky part!
        Num = 0.0;
        LINES = l;
        COLS = M-l+1;
        K.resize (LINES, COLS);
        for (i=0 ; i<M+1 ; i++) { // we "jump" above index j
          if (i<j) possible_indexes_for_x (i) = i;
          else if (i>j) possible_indexes_for_x (i-1) = i;
	}
        for (line=0 ; line<LINES ; line++) { 
          for (col=line ; col<COLS+line ; col++) K (line, col-line) = possible_indexes_for_x (col);
	}

        for (col=0 ; col<COLS ; col++) {
          k1_k2__kl.resize(LINES); // a vector of column indexes, one per line of K
          k1_k2__kl = col;
          while (k1_k2__kl (LINES - 2)>-1) {
            double prod_tmp = 1.0;
            for (line=0 ; line<LINES ; line++)  prod_tmp *= x (K (line, k1_k2__kl (line)));
            Num += prod_tmp;

            k1_k2__kl (0) = k1_k2__kl (0) - 1;
	    if ( (k1_k2__kl(0)==-1) && (LINES>2) ) {
              ind = count (k1_k2__kl (Range (0, LINES-3)) <= 0) - 1;
              k1_k2__kl (ind+1) = k1_k2__kl (ind+1) - 1;
              k1_k2__kl (Range (0, ind)) = k1_k2__kl (ind+1);
	    }

	  } // end while

	} // end for

      } // end else

      Num *= Sign;      

      // denominator
      for (i=0 ; i<M+1 ; i++) { // we "jump" above x (j)
        if (i<j) Den_tmp (i) = x (j) - x (i);
	else if (i>j) Den_tmp (i-1) = x (j) - x (i);
      }
      Den = product (Den_tmp);

      S_x(j, M-l) = Num/Den;

    } // for l
  } // for j
}

// test function
/*int main (void) {

  Array<double, 1> x (8);
  x = 1, 2, 3, 4, 7, 9, 13, 8;
  Array<double, 2> S_x;
  AIM_S (S_x, x);
  cout << S_x << endl;
  return 0; 
}
*/
