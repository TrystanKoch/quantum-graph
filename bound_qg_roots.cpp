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

// 
//
//
//
//

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





const double X_MIN = +0.1;
const double X_MAX = +10000;
const double Y_MIN = -0.01;
const double Y_MAX = +0.02; 


// Number of discrete points along a side
// 2**10 + 1 -> 10 halvings possible
const unsigned int NUM_HALVINGS_X = 22;
const unsigned int NUM_HALVINGS_Y = 8;


int main()
{

  // Read the input file. Pass the header along to the output files and
  // then send the data to a vector for further processing.
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

  // Initialize the quantum graph.
  QuantumGraph QG = QuantumGraph(header);
  
  // Instantiate the vectors
  std::clog << std::endl;
  std::clog << "Making Initial Boundary..." << std::endl;
  BoundingBox initialBoundary;
  initialBoundary = initializeBoundingBox(X_MIN, 
                                          X_MAX,
                                          Y_MIN,
                                          Y_MAX,
                                          NUM_HALVINGS_X,
                                          NUM_HALVINGS_Y, QG);

  // Now all of the boundary arguments are held in their proper vectors.
  // From now on, we can sum these vectors. The sum should be the number
  // of zeros inside the boundary.
  std::clog << std::endl;
  std::clog << "Initial Boundary Made." << std::endl;
  std::clog << "Total Winding Number: " << std::flush;
  std::clog << windingNumber(initialBoundary) << std::endl;
  std::clog << std::endl;



  //////////////////////////////////////////////////////////////////////
  /// Test the Recursive RootBounding
  
  std::cout << QG << std::endl;
  // We want to test the recursive root bounding algorithm. So 
  std::clog << "Recursively Finding Complex Root Bounds" << std::endl;
  unsigned long long int total 
      = recursiveRootBounding(initialBoundary, QG);

  std::clog << std::endl;
  std::clog << "Total Eigenvalues Bounded: " << total << std::endl;
  std::clog << std::endl;

  return 0;

}



