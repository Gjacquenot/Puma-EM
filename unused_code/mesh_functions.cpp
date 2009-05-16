/***************************************************************************
 * mesh_functions.cpp  some mesh functions
 *            
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


void edges_numbers_triangles_indexes(Array<int, 1>& indexes_of_triangles, const Array<int, 1>& list_of_edges_numbers, const Array<int, 2>& edges_numbers_triangles) {
  /* This function returns a 1-D array of the indexes of the triangles corresponding
   * to a 1-D array of edges_numbers. This function is important for creating lists of triangles
   * that will participate to the MoM or to the FMM for a given set of edges.
   */
  int i, j, tr_index;
  Range all = Range::all();
  const int N_edges = list_of_edges_numbers.size(), N_adj_triangles = edges_numbers_triangles.columns();

  // we first find the maximum and minimum of the triangles indexes
  int MIN_INDEX_TR = edges_numbers_triangles(0, 0), MAX_INDEX_TR = edges_numbers_triangles(0, 0);
  for (i=0 ; i<N_edges ; i++) {
    for (j=0 ; j<N_adj_triangles ; j++) {
      tr_index = edges_numbers_triangles(list_of_edges_numbers(i), j);
      if (tr_index == -1) break;
      else {
        if (tr_index<MIN_INDEX_TR) MIN_INDEX_TR = tr_index;
        if (tr_index>MAX_INDEX_TR) MAX_INDEX_TR = tr_index;
      }
    }
  }

  // we then fill an array containing one occurence of the indexes of the triangles
  Array<int, 1> indexes_of_triangles_tmp1(MAX_INDEX_TR - MIN_INDEX_TR + 1);
  indexes_of_triangles_tmp1 = -1;
  for (i=0 ; i<N_edges ; i++) {
    for (j=0 ; j<N_adj_triangles ; j++) {
      tr_index = edges_numbers_triangles(list_of_edges_numbers(i), j);
      if (tr_index == -1) break;
      else indexes_of_triangles_tmp1(tr_index-MIN_INDEX_TR) = tr_index;
    }
  }

  // we then compute the size of indexes_of_triangles and then fill it 
  int N_triangles = count( indexes_of_triangles_tmp1>-1 );
  indexes_of_triangles.resize(N_triangles);
  j = 0;
  for (i=0 ; i<indexes_of_triangles_tmp1.size() ; i++) {
    if (indexes_of_triangles_tmp1(i) > -1) {
      indexes_of_triangles(j) = indexes_of_triangles_tmp1(i);
      j++;
    }
  }
}


void computation_edges_numbers_local_edges_numbers(Array<int, 1>& edges_numbers_local_edges_numbers, const Array<int, 1>& list_of_edges_numbers, const Array<int, 2>& triangles_edges_numbers, const Array<int, 1>& indexes_of_triangles) {

  /* We must give "local" numbers to the edges, such that:
   * 
   *       local_edge_number = edges_numbers_local_edges_numbers[edge_number]
   * 
   * This "edges_numbers_local_edges_numbers" will be used in FMM and
   * MoM. This indirection allows for the creation of local matrices with smaller
   * dimensions than what should be expected by using directly the
   * global edges numbers.
   *    
   * The size of "edges_numbers_local_edges_numbers" is equal to the maximum edge number
   * encountered in the "triangles_edges_numbers[indexes_triangles]" 2-D array
   */
  Range all = Range::all();
  int i, MAX_EDGE_NUMBER = 0;
  const int N_edges = list_of_edges_numbers.size(), N_triangles = indexes_of_triangles.size();
  for (i=0 ; i<N_triangles ; i++) {
    int local_max = max( triangles_edges_numbers(indexes_of_triangles(i), all) );
    if (local_max > MAX_EDGE_NUMBER) MAX_EDGE_NUMBER = local_max;
  }

  edges_numbers_local_edges_numbers.resize(MAX_EDGE_NUMBER + 1);
  edges_numbers_local_edges_numbers = -1;

  // creation of the local numerotation; goes from 0 to N_edges-1
  for (i=0 ; i<N_edges ; i++) edges_numbers_local_edges_numbers(list_of_edges_numbers(i)) = i;
}


