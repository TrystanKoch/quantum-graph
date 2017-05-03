/*----------------------------------------------------------------------
   openquantumgraph.h
     - Class definition and Method prototypes for a class representing 
       the structure of a quantum graph.
     - Class Definitions and Function prototypes for classes modeling
       nodes, bonds, and undirected bonds in the quantum graph's 
       structure.

   Copyright (C) 2016  Trystan Koch

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------*/

#include "quantumgraph.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

// Use this code with the flag -DLAPACK_COMPLEX_CUSTOM. This way, the
// LAPACKE code can take the gsl_complex types. Because of the rather
// universal way that complex types are defined, this is desirable for
// removing much complexity from the code.
#define lapack_complex_double gsl_complex
#define lapack_complex_float gsl_complex_float

#include <lapacke.h>


OpenQuantumGraph::OpenQuantumGraph()
{
  
}





































