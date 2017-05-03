/*----------------------------------------------------------------------
   bound_qg_roots.cpp 
     - Using the argument principle, recursively bounds complex roots of
       a Quantum Graph. Outputs these bounds.
   Copyright (C) 2015  Trystan Koch

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

// In the spectral theory of Quantum Graphs, we often wish to find the
// eigenvalues of the graph. These eigenvalues are roots of the graph's
// characteristic function. For a lossless graph the characteristic function
// can be turned into a real function, however this involves additional
// operations in general, which is not desired. Furthermore, for a lossy graph
// this function has complex roots.
//
// Because complex roots exist off of the real number line, we cannot locate one
// as easily as we could in a real function with only real arguments. While
// Newton's method does work for complex functions, the basin of attraction for
// any given root is fractal for most. Since our problem has theoretically an
// infinite number of roots, with apparently random displacements from each
// other, we could never be sure that we would find all of the roots in an
// interval using only Newton's method. In addition, we would have a great deal
// of redundant roots, found from multiple starting points; even if we deleted 
// the useless ones, roots with higher multiplicity would appear the same as
// those of multiplicity one.
//
// So, in order to find complex roots in general, we must bound the roots first
// and then, with sufficiently small bounding boxes, we may initiate a search
// using Newton's method.
//
// This program uses Cauchy's Argument Principle (as used in Nyquist's stability
// theorem) to count the number of roots within a chosen box in the 
// complex plane. There is a more detailed discussion of this process in
//    quantumgraphrootfinding.h
//    quantumgraphrootfinding.cpp
// which define methods and classes for that calculation.
//
// This program takes in a quantum graph, creates the box to search for roots,
// determines how many roots are in the box, then calls the recursive function
// to isolate each root in its own box. Finally, it logs how many roots have 
// been bounded in this way as a check, to compare with the first number found.
// 
// The recursive root-finding function outputs all of the root bounding boxes 
// to standard output. Be aware of this before running!


#include "quantumgraph.h"
#include "quantumgraphrootfinding.h"

#include <cmath>
#include <limits>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <iostream>




const double X_MIN = +5;
const double X_MAX = +20;

const double DX = 1E-3;




int main()
{

  // Reads the input file. Passes the header along to the output files and
  // then sends the data to a vector for further processing.
  std::string line;
  std::vector<std::string> header;
  while (std::getline(std::cin, line))
  {
    // Check if there is a header.
    if (line[0]=='#')
    {
      header.push_back(line);
      continue;
    }
  }
  

  // Initializes the quantum graph from the input file.
  QuantumGraph QG = QuantumGraph(header);
  
  std::cout << QG << std::endl;
  
  gsl_complex z = gsl_complex_rect(X_MIN, 0.0);
  gsl_complex determinant;
  while (GSL_REAL(z) < X_MAX)
  {
    std::cout << GSL_REAL(z) << '\t';
    
    determinant = QG.RealCharacteristic(z);
    
    std::cout << GSL_REAL(determinant) << '\t';
    std::cout << GSL_IMAG(determinant) << std::endl;
    
    z = gsl_complex_add_real(z, DX);
  }
  

  return 0;

}



