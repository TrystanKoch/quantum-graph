/*----------------------------------------------------------------------
   ComplexQGRootFinder.cpp - Source for programs to find zeros of the
                           Quantum Graph characteristic equation in the
                           complex plane, using the argument principle
                           and Newton's Method.
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


#include "quantumgraph.h"
#include "quantumgraphrootfinding.h"

#include <limits>
#include <iostream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_qrng.h>





BoundingBox initializeBoundingBox(double realMin, double realMax, 
                                    double imagMin, double imagMax,
                                    unsigned int realHalvings, 
                                    unsigned int imagHalvings,
                                    const QuantumGraph &QG)
{
  BoundingBox BB;

  BB.realHalvings = realHalvings;
  BB.imagHalvings = imagHalvings;

  BB.realMin = realMin;
  BB.realMax = realMax;
  BB.imagMin = imagMin;
  BB.imagMax = imagMax;

  BB.bottom = makeArgSumVector(Bottom, imagMin, realMin, realMax,
                               realHalvings, QG);
  BB.right = makeArgSumVector(Right, realMax, imagMin, imagMax,
                              imagHalvings,QG);
  BB.top = makeArgSumVector(Top, imagMax, realMax, realMin,
                            realHalvings,QG);
  BB.left = makeArgSumVector(Left, realMin, imagMax, imagMin,
                             imagHalvings, QG);

  return BB;
}


//
//
//
//
//
std::ostream& operator<<(std::ostream& os, const BoundingBox BB)
{
  os.precision(std::numeric_limits<double>::max_digits10);
  os << std::scientific;
  os << std::showpos;

  os << (BB.realMin) << "    " << (BB.realMax) << "    " 
     << (BB.imagMin) << "    " << (BB.imagMax) << std::flush;

  return os;
}


// Given the argument sum vectors in a bounding box, tells us the
// winding number.
unsigned long long int windingNumber(const BoundingBox& BB)
{
  return (unsigned long long) ((( BB.bottom.back() - BB.bottom.front()
                                  + BB.right.back() - BB.right.front()
                                  + BB.top.back() - BB.top.front()
                                  + BB.left.back() - BB.left.front() ) 
                                / (2*M_PI)) + 0.5);
}


// Creates a vector along side "side" of the square boundary of the
// region
std::vector<double> makeArgSumVector(SideType side, double sideConstant,
                                     double sideStart, double sideEnd,
                                     unsigned int halvingsRemaining,
                                     const QuantumGraph &QG)
{
  // The difference in the argument at the begining from the argument at
  // The beginning is always zero.
  std::vector<double> argSumVector = {0};
  unsigned long long int numPoints = (1 << halvingsRemaining) + 1;
  double spacing = (sideEnd - sideStart) / (numPoints - 1);
  double argDiffSum = 0;

  gsl_complex z;
  gsl_complex fOld;
  gsl_complex fNew;

  if (side == Bottom || side == Top)
  {
    GSL_SET_COMPLEX(&z, sideStart, sideConstant);
    fOld = QG.Characteristic(z);
    for (unsigned long long int i=1; i<numPoints; i++)
    {
      GSL_SET_REAL(&z , sideStart + i * spacing);
      fNew = QG.Characteristic(z);
      argDiffSum += gsl_complex_arg(gsl_complex_div(fNew, fOld));
      argSumVector.push_back(argDiffSum);
      fOld = fNew;
    }
    return argSumVector;
  }
  else if (side == Right || side == Left)
  {
    GSL_SET_COMPLEX(&z, sideConstant, sideStart);
    fOld = QG.Characteristic(z);
    for (unsigned long long int i=1; i<numPoints; i++)
    {
      GSL_SET_IMAG(&z , sideStart + i * spacing);
      fNew = QG.Characteristic(z);
      argDiffSum += gsl_complex_arg(gsl_complex_div(fNew, fOld));
      argSumVector.push_back(argDiffSum);
      fOld = fNew;
    }
    return argSumVector;
  }
  else
  {
    return {};
  }
}


// When it reaches a box with winding number 1, prints the bounds to
// 
//
//
//
unsigned long long int 
recursiveRootBounding(BoundingBox& BB, const QuantumGraph &QG)
{
  unsigned long long int winding = windingNumber(BB);

  if (winding > 1)
  {
    if ((BB.realHalvings>1)&&(BB.imagHalvings>1))
    {
      if (((BB.realMax)-(BB.realMin)) > ((BB.imagMax)-(BB.imagMin)))
      {
        BoundingBox rightBB;
        splitBoxLeftRight(BB, rightBB, QG);
        return recursiveRootBounding(BB, QG) 
               + recursiveRootBounding(rightBB, QG);
      }
      else
      {
        BoundingBox bottomBB;
        splitBoxTopBottom(BB, bottomBB, QG);
        return recursiveRootBounding(BB, QG)
               + recursiveRootBounding(bottomBB, QG);
      }
    }
    else if (BB.realHalvings>1)
    {
      BoundingBox rightBB;
      splitBoxLeftRight(BB, rightBB, QG);
      return recursiveRootBounding(BB, QG) 
             + recursiveRootBounding(rightBB, QG);
    }
    else if (BB.imagHalvings>1)
    {
        BoundingBox bottomBB;
        splitBoxTopBottom(BB, bottomBB, QG);
        return recursiveRootBounding(BB, QG)
               + recursiveRootBounding(bottomBB, QG);
    }
    else
    {
      std::cout << BB << "    " << winding << std::endl;
      std::clog << BB << "    " << winding << std::endl;
    }
  }
  else if (winding == 1)
  {
    std::cout << BB << std::endl;
  }

  return winding;
}



//
//
//
//
//
void splitBoxLeftRight(BoundingBox& leftBB, BoundingBox& rightBB,
                       const QuantumGraph &QG)
{

  // The next set of boxes has been halved in the real dimension
  leftBB.realHalvings -= 1;
  rightBB.realHalvings = leftBB.realHalvings;
  rightBB.imagHalvings = leftBB.imagHalvings;

  // The "right" vector is the right vector of the right hand box
  leftBB.right.swap(rightBB.right);

  // Nothing about the imaginary limits have been changed
  rightBB.imagMax = leftBB.imagMax;
  rightBB.imagMin = leftBB.imagMin;
  rightBB.realMax = leftBB.realMax;

  // A new line at the center of the box must be in the new boundaries.
  // So the minimum and maximum of the real part change.
  double center = (leftBB.realMax + leftBB.realMin) / 2.0;
  leftBB.realMax = center;
  rightBB.realMin = center;

  // The latter half of the top vector goes to the left box while
  // the first half goes to the right box.
  leftBB.top.swap(rightBB.top);
  leftBB.top.assign(
      std::make_move_iterator(rightBB.top.begin() 
                              + (1<<rightBB.realHalvings)),
      std::make_move_iterator(rightBB.top.end()));

  rightBB.top.erase(
      rightBB.top.begin() + ((1<<rightBB.realHalvings) + 1) ,
      rightBB.top.end());


  // The latter half of the bottom vector goes to the right box while
  // the first half goes to the left box.
  rightBB.bottom.assign (
      std::make_move_iterator(leftBB.bottom.begin()
                              + (1<<leftBB.realHalvings)),
      std::make_move_iterator(leftBB.bottom.end()));

  leftBB.bottom.erase(
      leftBB.bottom.begin() + ((1<<leftBB.realHalvings) + 1),
      leftBB.bottom.end());

  // Finally, make the center boundary line for the two sides
  leftBB.right = makeArgSumVector(Right, center, 
                                  leftBB.imagMin, leftBB.imagMax, 
                                  leftBB.imagHalvings, QG);
  rightBB.left = makeArgSumVector(Left, center, 
                                  leftBB.imagMax, leftBB.imagMin, 
                                  rightBB.imagHalvings, QG);
}



//
//
//
//
//
void splitBoxTopBottom(BoundingBox &topBB, BoundingBox &bottomBB, 
                       const QuantumGraph &QG)
{

  // The next set of boxes has been halved in the real dimension
  topBB.imagHalvings -= 1;
  bottomBB.realHalvings = topBB.realHalvings;
  bottomBB.imagHalvings = topBB.imagHalvings;

  // The "bottom" vector is the bottom vector of the bottom box
  topBB.bottom.swap(bottomBB.bottom);

  // Nothing about the real limits have been changed
  bottomBB.realMax = topBB.realMax;
  bottomBB.realMin = topBB.realMin;
  bottomBB.imagMin = topBB.imagMin;

  // A new line at the center of the box must be in the new boundaries.
  // So the minimum and maximum of the real part change.
  double center = (topBB.imagMax + topBB.imagMin) / 2.0;
  topBB.imagMin = center;
  bottomBB.imagMax = center;

  // The latter half of the right vector goes to the top box while
  // the first half goes to the bottom box.
  topBB.right.swap(bottomBB.right);
  topBB.right.assign(
      std::make_move_iterator(bottomBB.right.begin() 
                              + (1<<bottomBB.imagHalvings)),
      std::make_move_iterator(bottomBB.right.end()));

  bottomBB.right.erase(
      bottomBB.right.begin() + ((1<<bottomBB.imagHalvings) + 1) ,
      bottomBB.right.end());


  // The latter half of the left vector goes to the bottom box while
  // the first half goes to the top box.
  bottomBB.left.assign (
      std::make_move_iterator(topBB.left.begin()
                              + (1<<topBB.imagHalvings)),
      std::make_move_iterator(topBB.left.end()));

  topBB.left.erase(
      topBB.left.begin() + ((1<<topBB.imagHalvings) + 1),
      topBB.left.end());

  // Finally, make the center boundary line for the two sides
  topBB.bottom = makeArgSumVector(Bottom, center, 
                                  topBB.realMin, topBB.realMax, 
                                  topBB.realHalvings, QG);
  bottomBB.top = makeArgSumVector(Top, center, 
                                  bottomBB.realMax, bottomBB.realMin, 
                                  topBB.realHalvings, QG);
}

//
//
//
//
gsl_complex newtonRootFinderQRNG(BoundingBox& BB, gsl_qrng * qrng,
                                   const QuantumGraph &QG)
{
  unsigned int restart = 0;
  gsl_complex zOld, zNew;
  double x[2];

  gsl_qrng_init(qrng);
  gsl_qrng_get(qrng, x);
  GSL_SET_COMPLEX(&zNew, BB.realMin*(1-x[0]) + x[0]*BB.realMax,
                         BB.imagMin*(1-x[1]) + x[1]*BB.imagMax);
  do
  {
    zOld = zNew;  
    zNew = QG.NewtonStep(zOld);
    if (outOfBoundary(zNew,BB))
    {
      gsl_qrng_get(qrng, x);
      GSL_SET_COMPLEX(&zNew, BB.realMin*(1-x[0]) + x[0]*BB.realMax,
                             BB.imagMin*(1-x[1]) + x[1]*BB.imagMax);
      restart += 1;
      /*std::clog.precision(17);
      std::clog << std::scientific << std::showpos;
      std::clog << std::endl;
      std::clog << restart <<std::endl;
      std::clog << BB.realMin << " "<<BB.imagMin <<std::endl;
      std::clog << x[0]<< " "<<x[1] <<std::endl;*/

    }
    else if (gsl_complex_abs(gsl_complex_sub(zNew,zOld))/gsl_complex_abs(zOld) <NEWTON_QG_EPS)
    {/*
 BB =  initializeBoundingBox(BB.realMin, BB.realMax, 
                              BB.imagMin, BB.imagMax,
                              NEWTON_INIT_HALVINGS,
                              NEWTON_INIT_HALVINGS,
                              QG);
  std::clog << "1 Root found in [" << BB << "]" << std::endl;
  std::clog << "   Winding Number: " << windingNumber(BB) << std::endl;//*/
      return zNew;
    }
  //std::clog << GSL_REAL(zOld)<< "  " <<GSL_IMAG(zOld) << std::endl;
  
  } while (restart <= NEWTON_QG_MAX_REPEATS);

  //New, 9/24
  BB =  initializeBoundingBox(BB.realMin, BB.realMax, 
                              BB.imagMin, BB.imagMax,
                              NEWTON_INIT_HALVINGS,
                              NEWTON_INIT_HALVINGS,
                              QG);
  std::clog << "No root found in [" << BB << "]" << std::endl;
  std::clog << "   Winding Number: " << windingNumber(BB) << std::endl;
  return GSL_COMPLEX_ZERO;
}

/*
gsl_complex newtonRootFinderRecursive(BoundingBox&,
                                      const QuantumGraph&)
{
  unsigned int restart = 0;
  gsl_complex zOld, zNew;

  GSL_SET_COMPLEX(&zNew, BB.realMin*(1-x[0]) + x[0]*BB.realMax,
                         BB.imagMin*(1-x[1]) + x[1]*BB.imagMax);
  do
  {
    zOld = zNew;  
    zNew = QG.newtonStep(zOld);
    if (outOfBoundary(zNew,BB))
    {
      GSL_SET_COMPLEX(&zNew, BB.realMin*(1-x[0]) + x[0]*BB.realMax,
                             BB.imagMin*(1-x[1]) + x[1]*BB.imagMax);
      restart += 1;
    }
    else if (gsl_complex_abs(gsl_complex_sub(zNew,zOld))/gsl_complex_abs(zOld) <NEWTON_QG_EPS)
    {
      return zNew;
    }

  } while (restart <= NEWTON_QG_MAX_REPEATS);

  std::clog << "No root found in [" << BB << "]" << std::endl;
  return GSL_COMPLEX_ZERO;
}
//*/

bool outOfBoundary(gsl_complex z, const BoundingBox& BB)
{
  if (GSL_REAL(z) < BB.realMin || GSL_REAL(z) > BB.realMax
      || GSL_IMAG(z) < BB.imagMin || GSL_IMAG(z) > BB.imagMax)
  {
    return true;
  }
  else
  {
    return false;
  }
}


