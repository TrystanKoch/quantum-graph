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



// These constants define the boundary of the box in the complex plane.
// X denotes the real line, Y denotes the imaginary line.
//
// Right now this is the only way to specify a particular box. This may 
// change in the future if I figure out an elegent way to pass these four
// numbers in a way that is understandable.
const double X_MIN = +1000;
const double X_MAX = +1100;
const double Y_MIN = -0.011;
const double Y_MAX = +0.01;


// Number of discrete points along a side. Because we are using the same
// argument information over and over again, we want a number of points that
// will divide evenly by two, multiple times.
//
// So instead of declaring the number of points at which this program should
// calculate the argument, these constants define maximum number of recursive
// steps that the user wants to take in the real and imaginary boundaries.
//
// That is, if I specify NUM_HALVINGS_X=10, then I will calculate the quantum
// graph's characterstic function (2**10 + 1) times along the boundary with a
// constant imaginary part. 
//
// Choose carefully! There is not quite a way to ensure that you will calculate
// enough points that the winding number is correct. It should be clear from
// Weyl's law whether the original winding number is approximately right, so the
// user can restart the calculation if it gives a very different number. 
// However, especially for graphs that have Time reversal symmetry, and 
// therefore only linear repulsion between eigenvalues, that two very close
// eigenvalues might go unnoticed. (Still working on way to fix that.)
const unsigned int NUM_HALVINGS_X = 15;
const unsigned int NUM_HALVINGS_Y = 8;




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
  
  // Instantiates the initial boundary. Quite a bit of calculation goes on here.
  // The argument of the characteristic function must be found at each point
  // along the boundary. This is the most efficient way to keep the argument
  // information, which is then passed down recursively later.
  std::clog << std::endl;
  std::clog << "Making Initial Boundary..." << std::endl;
  BoundingBox initialBoundary;
  initialBoundary = initializeBoundingBox(X_MIN, 
                                          X_MAX,
                                          Y_MIN,
                                          Y_MAX,
                                          NUM_HALVINGS_X,
                                          NUM_HALVINGS_Y, QG);

  // Now all of the boundary arguments are held in their proper vectors within
  // a BoundingBox object. From now on, we can sum these vectors.
  // quantumgraphrootfinding.cpp contains the functions necessary to do that
  // work. The sum should be the number of zeros inside the boundary. 
  std::clog << std::endl;
  std::clog << "Initial Boundary Made." << std::endl;
  std::clog << "Total Winding Number: " << std::flush;
  std::clog << windingNumber(initialBoundary) << std::endl;
  std::clog << std::endl;



  //////////////////////////////////////////////////////////////////////
  /// Recursive RootBounding
  
  // Write the Quantum Graph information as a header to the output file. The
  // output must be redirected from standard output!
  std::cout << QG << std::endl;
  
  // Now Comes the recursive step. Depending only on the number of times you
  // have chosen to recursively halve the original bounding box, this function
  // calculates the winding number of the boundary, cuts the box in two,
  // calculates the arguments along the bisecting line, and passes all of the
  // information to the next recursive step (one new box at a time). It
  // continues until either it finds an empty box, a box that contains one and
  // only one root, or the recursive steps run out. 
  //
  // All boxes containing one root are sent to standard output. The function
  // returns the number of bounding boxes it found.
  std::clog << "Recursively Finding Complex Root Bounds" << std::endl;
  unsigned long long int total 
      = recursiveRootBounding(initialBoundary, QG);
        
  // Gives the user the number of bounding boxes found, as a comparison to the
  // winding number above.
  std::clog << std::endl;
  std::clog << "Total Eigenvalues Bounded: " << total << std::endl;
  std::clog << std::endl;

  return 0;

}



