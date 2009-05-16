/***************************************************************************
 * create_mesh.cpp  function that creates the mesh needed for the computations
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
#include <fstream>
#include <blitz/array.h>

using namespace blitz;

#include "mesh.h"

void create_mesh (mesh & MESH, const char * nodes_name, const char * triangles_name, const char * edges_name) {

  int i, j;
  Array<double, 2> nodes_tmp, triangles_tmp, edges_tmp;

  // reading nodes.txt into nodes_tmp
  {
    // char * filename = "nodes.txt";
    ifstream ifs (nodes_name);
    if (! ifs.is_open()) { 
      cout << "create_mesh : Error opening file : " << nodes_name << endl; 
      exit (1); 
    }
    ifs >> nodes_tmp;
    ifs.close();
  }

  // reading triangles.txt into triangles_tmp
  {
    // char * filename = "triangles.txt";
    ifstream ifs (triangles_name);
    if (! ifs.is_open()) { 
      cout << "create_mesh : Error opening file : " << triangles_name << endl; 
      exit (1); 
    }
    ifs >> triangles_tmp;
    ifs.close();
  }

  // reading edges.txt into edges_tmp
  {
    // char * filename = "edges.txt";
    ifstream ifs (edges_name);
    if (! ifs.is_open()) { 
      cout << "create_mesh : Error opening file : " << edges_name << endl; 
      exit (1); 
    }
    ifs >> edges_tmp;
    ifs.close();
  }

  // creation of MESH.edges
  MESH.edges.resize(edges_tmp.rows(), edges_tmp.columns()-1);
  MESH.edges = cast<int> (edges_tmp(Range::all(), Range(0, 5)));

  // creation of MESH.edge_length
  MESH.edge_length.resize(edges_tmp.rows());
  MESH.edge_length = edges_tmp(Range::all(), 6);

  // MESH.z_min, MESH.z_max, MESH.rho_max
  MESH.x_min = min (nodes_tmp (Range::all(), 1)), MESH.x_max = max (nodes_tmp (Range::all(), 1));
  MESH.y_min = min (nodes_tmp (Range::all(), 2)), MESH.y_max = max (nodes_tmp (Range::all(), 2));
  MESH.z_min = min (nodes_tmp (Range::all(), 3)), MESH.z_max = max (nodes_tmp (Range::all(), 3));
  MESH.rho_min = 0.0;
  MESH.rho_max = sqrt (pow2 (MESH.x_max-MESH.x_min) + pow2 (MESH.y_max-MESH.y_min));

  // creation of MESH.nodes
  // the non-attributed lines = -1
  int node_number_max = (int) max(nodes_tmp(Range::all(), 0));
  int N_nodes = nodes_tmp.rows(), node_number;
  MESH.nodes.resize(node_number_max+1, 4);
  MESH.nodes = 0.0;
  MESH.nodes(Range::all(), 0) = -1.0;
  for (j=0 ; j<N_nodes ; j++) {
    node_number = (int) (nodes_tmp(j, 0));
    MESH.nodes(node_number, Range::all()) = nodes_tmp(j, Range::all());
  }

  // creation of MESH.tr_numbers
  int N_triangles = triangles_tmp.rows(), index;
  MESH.tr_numbers.resize(N_triangles);
  MESH.tr_numbers = cast<int> (triangles_tmp(Range::all(), 0));

  // creation of MESH.triangles
  MESH.triangles.resize(N_triangles, 3);
  MESH.triangles = cast<int> (triangles_tmp(Range::all(), Range(1, 3)));

  // creation of MESH.triangles_E_H_J_M
  MESH.triangles_E_H_J_M.resize(N_triangles, 4);
  MESH.triangles_E_H_J_M = cast<int> (triangles_tmp(Range::all(), Range(4, 7)));

  // creation of MESH.ind_edges_of_tr_i
  int triangle_number, col;
  MESH.ind_edges_of_tr_i.resize(N_triangles, 4);
  MESH.ind_edges_of_tr_i = -1;
  for (j=0 ; j<N_triangles ; j++) {
    triangle_number = MESH.tr_numbers (j);
    col = 0;
    for (index = 0 ; index < MESH.edges.rows() ; index++) {
      if (MESH.edges (index, 0) == triangle_number) {
        MESH.ind_edges_of_tr_i (j, col) = index;
        col++;
      }
    }
  }
  //cout << MESH.ind_edges_of_tr_i << endl;
  //cout << MESH.edges << endl;

}

