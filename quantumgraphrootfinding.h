/*----------------------------------------------------------------------
   quantumgraphrootfinding.h 
     - Headers for programs to find zeros of the Quantum Graph 
       characteristic equation in the complex plane, using the argument
       principle and Newton's Method.

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


#ifndef COMPLEXQGROOTFINDER_H_
#define COMPLEXQGROOTFINDER_H_


#include "quantumgraph.h"

#include <iostream>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_qrng.h>




// We assume Newton's root-finding method converged if another iteration
// does not change the value of the output. This is the absolute value
// of the difference in the complex plane that the program considers to
// be "not change".
const double NEWTON_QG_EPS = 5e-16;





// We may have found a bounding box which is much larger than the basin
// of convergence around the root we wish to find. Though the qrng
// method attempts to sample the bounding box for this basin, it's not
// useful to spend too much effort looking for one root. If this many
// initial conditions do not end in Newton's method converging, then
// we output an error message and move on.
const unsigned int NEWTON_QG_MAX_REPEATS = 8000;



//
//
//
//
//
const unsigned int NEWTON_INIT_HALVINGS = 15;





// This enumeration makes it clear which side of the boundary we are
// creating at a given time. This way, I do not need to write four
// different vector creation functions for the vectors in BoundingBox.
//
// Not necessary, just makes the makeArgSumVector function slightly 
// easier to understand.
enum SideType {Bottom, Right, Top, Left};





// Represents a rectangular region in the complex plane. The real and
// imaginary minima and maxima are stored.
//
// Because the program I have written relies on recursively splitting
// the bounded region into halves, the struct does not contain the
// number of sampled points along the boundaries, but instead contains
// the number of times we can cut the box in half without evaluating any
// additional points along line segments we already used. That is, the
// number of points at which our program samples the function f(z) will
// be 2^N + 1, where N is "realHalvings" or "imagHalvings" depending on
// whether the line goes parallel or perpendicularly (respectively) to
// the real axis.
//
// Finally, each side of the boundary has a vector. Because we only care
// about the accumulated change in the argument of the function as we
// trace out a path in the complex plane around the boundary, we store
// the sum of the argument along the side. The total additional argument
// from moving from the beginning to the end of the boundary line is the
// last element of the vector minus the first element. This is still
// true if we cut the vector into two pieces.
//
// I chose not to make this a class for memory reasons. Because there
// will be many copies of this data passed recursively, I wanted to keep
// functions from being passed along as well.
struct BoundingBox
{
  unsigned int realHalvings;
  double realMin;
  double realMax;

  unsigned int imagHalvings;
  double imagMin;
  double imagMax;

  std::vector<double> bottom = {};
  std::vector<double> right = {};
  std::vector<double> top = {};
  std::vector<double> left = {};
};





// Creates a new BoundingBox, which is the rectangle in the complex
// plane that stretches from realMin to realMax in the real direction
// and imagMin to imagMax in the imaginary direction.
//
// The realHalvings and imagHalvings are related to the number of points
// that are sampled in the real and imaginary directions, respectively,
// along the boundary. The number of points is 2^N + 1, ensuring that we
// can always cleanly divide the boundary into two equally sized halves. 
//
// Because the BoundingBox struct's argument-sum vectors depend on the
// Quantum Graph I pass to it, this is also an argument. The function
// takes a constant reference to a quantum graph because it will not
// change the graph, only use it to find the arguments.
BoundingBox initializeBoundingBox(double realMin, double realMax, 
                                    double imagMin, double imagMax,
                                    unsigned int realHalvings, 
                                    unsigned int imagHalvings,
                                    const QuantumGraph &QG);




// Prints the real and imaginary bounds of the bounding box to a line
// in some ostream object (most likely std::cout). It prints 
//   realMin  realMax  imagMin  imagMax
// so that we can break the root-bounding and root-polishing to two
// different steps, seperated by a file containing the bounds for the
// complex roots.
//
// Not necessary, but cleans up the recursive function below.
std::ostream& operator<<(std::ostream&, const BoundingBox);





// The next two functions do basically the same thing. The first splits
// a bounding box into a left half and a right half, while the second
// splits a bounding box into a top and a bottom half.
//
// THE FIRST ARGUMENT IS BOTH AN OUTPUT AND AN INPUT. This is to save
// memory and time. It takes in a reference to the original bounding box
// and, after the function has run, the left (top) bounding box is now
// stored where the original was. We can do this because once the box is
// passed down the recursive chain, we never need it again elsewhere. We
// stop caring about the original box because we already determined it
// had more than one zero inside it and we are only interested in
// boundaries with only one zero inside.
//
//
//
// The right (bottom) bounding box must be initialized before splitting.
// A fair bit of memory swapping happens in this function, since it is
// the function that, after the initial box creation, is called over and
// over again in the recursion.
void splitBoxLeftRight(BoundingBox&, BoundingBox&,
                       const QuantumGraph&);
void splitBoxTopBottom(BoundingBox&, BoundingBox&, 
                       const QuantumGraph&);




// After initialization, a bounding box owns four vectors full of sums
// of the change in the argument as we go around the boundary from 
// sampled point to sampled point. This function takes the first and
// last elements from each of the four sides' argument-sum vectors and
// calculates the winding number around the boundary.
//
// While the vector stores the summed argument, the function calculates
// the winding number by dividing this by 2*pi. It's easier to do that
// division here than it is to do when we store the vector.
unsigned long long int windingNumber(const BoundingBox&);





//
//
//
//
//
std::vector<double> 
makeArgSumVector(SideType, double, double, double, unsigned int,
                 const QuantumGraph&);




//
//
//
//
//
unsigned long long int 
recursiveRootBounding(BoundingBox&, const QuantumGraph&);





// Once we've bounded the roots, we still need to find them. We can use
// Newton's Method to find them. Because Newton's method might converge
// to a different root, or travel off to infinity, if the root-finder
// leaves the bounding box on an iteration, we would like to choose
// another initial guess to use Newton's Method on. Because there will
// be eaxctly one root in each bounded box, we know that if the method
// converges within the box, we have found the correct root.
//
// This root-finder operates on the assumption that the best way to find
// an initial input to Newton's Method is to sample the bounding box
// quasi randomly. We assume that a calling function has set up a
// quasi-random number generator (I tend to use Sobol's algorithm)
// which it then initializes and calls until it finds an input for which
// Newton's Method converges within the box.
//
// 
//
// We assume that the zero we're trying to find has a basin of
// convergence that is comparable to the size of the bounding box. If
// this is not true, it might be possible for this function to 
// consistantly miss the basin and not return an answer. If it does this
// enough (I use 250 points, though it doesn't get nearly this large in
// practice) it cuts off and returns GSL_COMPLEX_ZERO, writing a warning
// to std::clog.
//
// In the case that a bounding box containing only one root is not small
// enough for this method to work, try cutting the size of the box using
// the argument principle.
gsl_complex newtonRootFinderQRNG(BoundingBox&, gsl_qrng* qrng,
                                 const QuantumGraph&);


// This Root-finder operates under the idea that we can continue to use
// the argument principle to reduce the size of the bounding box around
// the root. If the Newton Root-Finder does not converge to a root 
// within the bounding box, we halve the bounding box, and try again.
// We can ignore the 
gsl_complex newtonRootFinderRecursive(BoundingBox&,
                                      const QuantumGraph&);


// Checks to see if a complex number is outside of a given boundary.
//
//
//
// Really, just a helper function for the Newton's Method, since the
// zero has to be inside the isolating bounding box.
bool outOfBoundary(gsl_complex, const BoundingBox&);





#endif







