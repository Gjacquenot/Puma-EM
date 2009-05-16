/*
** mesh_functions.h
** Login : <denbosch@localhost.localdomain>
** Started on  Thu Dec 15 15:58:34 2005 Idesbald van den Bosch
** $Id$
** 
** Copyright (C) 2005 Idesbald van den Bosch
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

using namespace blitz;

void edges_numbers_triangles_indexes(Array<int, 1>& indexes_of_triangles, const Array<int, 1>& list_of_edges_numbers, const Array<int, 2>& edges_numbers_triangles);

void computation_edges_numbers_local_edges_numbers(Array<int, 1>& edges_numbers_local_edges_numbers, const Array<int, 1>& list_of_edges_numbers, const Array<int, 2>& triangles_edges_numbers, const Array<int, 1>& indexes_of_triangles);
