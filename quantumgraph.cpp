/*----------------------------------------------------------------------
   quantumgraph.cpp - Methods for a class representing a quantum graph
                      that includes loss.
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

// See header file for information regarding the member functions. This
// file contains all the source code for methods used by QuantumGraph
// objects.

#include "quantumgraph.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <limits>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

// Use this code with the flag -DLAPACK_COMPLEX_CUSTOM. This way, the
// LAPACKE code can take the gsl_complex types. Because of the rather
// universal way that complex types are defined, this is desirable for
// removing much complexity from the code.
#define lapack_complex_double gsl_complex
#define lapack_complex_float gsl_complex_float

#include <lapacke.h>

// The maximum value we are willing to accept from a singular value
// before we give a warning about matrix not being singular.
const double SVD_EPSILON = 1e-10;

////////////////////////////////////////////////////////////////////////
// Constructor Methods


QuantumGraph::QuantumGraph()
{
  numBonds = 1;

  allocGraphMemory(numBonds);

  makeInternals();
}

// Constructor from GSL types
QuantumGraph::QuantumGraph(const gsl_vector_complex* LL,
                           const gsl_matrix_complex* SS)
{
  numBonds = SS->size1;

  if (LL->size != numBonds)
  {
    std::clog << "Length Vector and S Matrix must be same dimension." 
              << std::endl;
  }
  else if (SS->size2 != numBonds)
  {
    std::clog << "S Matrix must be square." 
              << std::endl;
  }

  allocGraphMemory(numBonds);

  gsl_matrix_complex_memcpy(SMatrix, SS);
  gsl_vector_complex_memcpy(LengthVector, LL);

  makeInternals();
}




// Constructor from vectors of doubles representing real and imaginary
// parts of LengthVector and SMatrix.
QuantumGraph::QuantumGraph(const std::vector<double>LLreal, 
                           const std::vector<double>LLimag,
                           const std::vector<std::vector<double>> SSreal,
                           const std::vector<std::vector<double>> SSimag)
{
  // Make sure we're not giving the object vectors and matrices of
  // different sizes.
  numBonds = LLreal.size();
  if (numBonds!=LLimag.size())
  {
    std::clog << "Imaginary part of Length Vector has incorrect size." 
              << std::endl;
  }
  if (!checkMatrixSize(SSreal))
  {
    std::clog << "Imaginary part of S Matrix has incorrect size." 
              << std::endl;
  }
  if (!checkMatrixSize(SSimag))
  {
    std::clog << "Real part of S Matrix has incorrect size." 
              << std::endl;
  }

  // Okay, allocate the memory for data members.
  allocGraphMemory(numBonds);
  
  // Make SMatrix
  for (unsigned int i=0; i<numBonds; i++)
  {
      for (unsigned int j=0; j<numBonds; j++)
      {
          gsl_complex Sij = gsl_complex_rect(SSreal[i][j],SSimag[i][j]);
          gsl_matrix_complex_set(SMatrix, i, j, Sij);
      }
  }

  // Make LengthVector
  for (unsigned int j=0; j<numBonds; j++)
  {
      gsl_complex Lj = gsl_complex_rect(LLreal[j], LLimag[j]);
      gsl_vector_complex_set(LengthVector, j, Lj);
  }

  makeInternals();
}


// Constructor that examines a file header.
//  Each program's output has a header that contains the details of the
//  quantum graph. This constructor is for reading a string containing
//  that header and initializing a Quantum Graph based on it.
//  See the header file for information about the formatting.
QuantumGraph::QuantumGraph(std::vector<std::string> header)
{
  std::vector<std::vector<double>> SMatrixReal;
  std::vector<std::vector<double>> SMatrixImag;
  std::vector<double> LengthsImag;
  std::vector<double> LengthsReal;

  numBonds = 0;

  for (unsigned int i=0; i<header.size(); i++)
  {
    std::string line = header[i];
    if (line == "#    Number of Bonds:")
    {
      line = header[i+1];
      line.erase(0,1);
      numBonds = std::stoi(line);
      LengthsImag.resize(numBonds);
      LengthsReal.resize(numBonds);
      SMatrixImag.resize(numBonds);
      SMatrixReal.resize(numBonds);
      for (unsigned int m=0; m<numBonds; m++)
      {
        SMatrixImag[m].resize(numBonds);
        SMatrixReal[m].resize(numBonds);
      }
    }
    else if (line == "#    S-Matrix (Real Part):")
    {
      for (unsigned int m=0; m<numBonds; m++)
      {
        line = header[i+1+m];
        line.erase(0,1);
        std::stringstream SMatrixRealStream(line);
        for (unsigned n=0; n<numBonds; n++)
        {
          SMatrixRealStream >> SMatrixReal[m][n];
        }
      }
    }
    else if (line == "#    S-Matrix (Imaginary Part):")
    {
      for (unsigned int m=0; m<numBonds; m++)
      {
        line = header[i+m+1];
        line.erase(0,1);
        std::stringstream SMatrixImagStream(line);

        for (unsigned n=0; n<numBonds; n++)
        {
          SMatrixImagStream >> SMatrixImag[m][n];
        }
      }
    }
    else if (line == "#    Lengths (Real Part):")
    {
      line = header[i+1];
      line.erase(0,1);
      std::stringstream LengthRealStream(line);
      for (unsigned int m=0; m<numBonds; m++)
      {
        LengthRealStream >> LengthsReal[m];
      }
    }
    else if (line == "#    Lengths (Imaginary Part):")
    {
      line = header[i+1];
      line.erase(0,1);
      std::stringstream LengthImagStream(line);
      for (unsigned int m=0; m<numBonds; m++)
      {
        LengthImagStream >> LengthsImag[m];
      }
    }
    else
    {
      continue;
    }
  }

  if (numBonds!=LengthsReal.size())
  {
    std::clog << "Real part of Length Vector has incorrect size." 
              << std::endl;
  }
  if (numBonds!=LengthsImag.size())
  {
    std::clog << "Imaginary part of Length Vector has incorrect size." 
              << std::endl;
  }
  if (!checkMatrixSize(SMatrixReal))
  {
    std::clog << "Imaginary part of S Matrix has incorrect size." 
              << std::endl;
  }
  if (!checkMatrixSize(SMatrixImag))
  {
    std::clog << "Real part of S Matrix has incorrect size." 
              << std::endl;
  }
  

  allocGraphMemory(numBonds);

  for (unsigned int i=0; i<numBonds; i++)
  {
    for (unsigned int j=0; j<numBonds; j++)
    {
      gsl_complex Sij = gsl_complex_rect(SMatrixReal[i][j],
                                         SMatrixImag[i][j]);
      gsl_matrix_complex_set(SMatrix, i, j, Sij);
    }
  }

  for (unsigned int i=0; i<numBonds; i++)
  {
    gsl_complex Li = gsl_complex_rect(LengthsReal[i],
                                      LengthsImag[i]);
    gsl_vector_complex_set(LengthVector, i, Li);
  }

  makeInternals();
}




// Copy constructor.
QuantumGraph::QuantumGraph(const QuantumGraph& QG)
{
  numBonds = QG.numBonds;
  allocGraphMemory(numBonds);

  gsl_matrix_complex_memcpy(SMatrix, QG.SMatrix);
  gsl_vector_complex_memcpy(LengthVector, QG.LengthVector);

  makeInternals();
}




// Need to get rid of allocated data when we destroy the object, because
// GSL uses pure C and needs to allocate memory.
QuantumGraph::~QuantumGraph()
{
  freeGraphMemory();
}




////////////////////////////////////////////////////////////////////////
// Internal Data initialization Methods


void QuantumGraph::makeInternals()
{
  makeNegImaginaryLengthVector();
  makePosImaginaryLengthVectorSum();
}




// So we don't have to do this on every computation.
void QuantumGraph::makeNegImaginaryLengthVector()
{
  gsl_vector_complex_memcpy(negImaginaryLengthVector, LengthVector);
  gsl_vector_complex_scale(negImaginaryLengthVector, 
                           gsl_complex_rect(0, -1) );
}




// So we don't have to do this on every computation.
void QuantumGraph::makePosImaginaryLengthVectorSum()
{
  gsl_complex LengthSum = GSL_COMPLEX_ZERO;
  for (unsigned int i=0; i<numBonds; i++)
  {
    gsl_complex Li = gsl_vector_complex_get(LengthVector, i);
    LengthSum = gsl_complex_add(LengthSum, Li);
  }
  posImaginaryLengthVectorSum 
    = gsl_complex_mul(LengthSum, gsl_complex_rect(0, 1));
}

// Check the size of a vector of vectors during initialization.
bool QuantumGraph::checkMatrixSize(std::vector<std::vector<double>> M)
{
  if (M.size() != numBonds)
  {
    return 0;
  }
  for (unsigned int i=0; i<numBonds; i++)
  {
    if (M[i].size() != numBonds)
    {
      return 0;
    }
  }
  return 1;
}




////////////////////////////////////////////////////////////////////////
// Memory Handling Methods


void QuantumGraph::allocGraphMemory(unsigned int N)
{
  SMatrix = gsl_matrix_complex_calloc(N,N);
  LengthVector = gsl_vector_complex_calloc(N);
  negImaginaryLengthVector = gsl_vector_complex_calloc(N);
}




void QuantumGraph::freeGraphMemory()
{
  gsl_matrix_complex_free(SMatrix);
  gsl_vector_complex_free(LengthVector);
  gsl_vector_complex_free(negImaginaryLengthVector);
}

////////////////////////////////////////////////////////////////////////
/// Methods for Derived Classes only

// Protected method for resetting the representation of the graph based
// on the implementation of a derived class.
void QuantumGraph::setGraph(const gsl_vector_complex* LL,
                            const gsl_matrix_complex* SS)
{
  freeGraphMemory();

  numBonds = SS->size1;

  if (LL->size != numBonds)
  {
    std::clog << "Length Vector and S Matrix must be same dimension." 
              << std::endl;
  }
  else if (SS->size2 != numBonds)
  {
    std::clog << "S Matrix must be square." 
              << std::endl;
  }

  allocGraphMemory(numBonds);

  gsl_matrix_complex_memcpy(SMatrix, SS);
  gsl_vector_complex_memcpy(LengthVector, LL);

  makeInternals();
}


////////////////////////////////////////////////////////////////////////
// Computation Methods


// While this function will have zeros, it will not be strictly real
// along the real axis. Only a complex zero-finding algorithm will work.
gsl_complex QuantumGraph::characteristic(const gsl_complex z) const
{
  gsl_matrix_complex* D = gsl_matrix_complex_calloc(numBonds, numBonds);
  gsl_matrix_complex_memcpy(D, SMatrix);

  for (unsigned int j=0; j<numBonds; j++)
  {
    gsl_complex negILj
      = gsl_vector_complex_get(negImaginaryLengthVector, j);
    gsl_complex Tj = gsl_complex_exp(gsl_complex_mul(negILj, z));
    gsl_matrix_complex_set(D, j, j, gsl_complex_negative(Tj));
  }

  int signum;
  gsl_permutation* p = gsl_permutation_alloc(numBonds);

  gsl_linalg_complex_LU_decomp(D, p, &signum);
  gsl_complex determinant = gsl_linalg_complex_LU_det(D, signum);

  gsl_permutation_free(p);
  gsl_matrix_complex_free(D);

  return determinant;
}





gsl_vector_complex* QuantumGraph::eigenvector(const gsl_complex z) const
{
  gsl_matrix_complex* D = gsl_matrix_complex_calloc(numBonds, numBonds);
  gsl_matrix_complex_memcpy(D, SMatrix);

  for (unsigned int j=0; j<numBonds; j++)
  {
    gsl_complex negILj
      = gsl_vector_complex_get(negImaginaryLengthVector, j);
    gsl_complex Tj = gsl_complex_exp(gsl_complex_mul(negILj, z));
    gsl_matrix_complex_set(D, j, j, gsl_complex_negative(Tj));
  }
  
  // The GSL library does not offer a complex singular value
  // decomposition routine. So we must rely on LAPACK, which has a
  // standard c routine interface called LAPACKE.
  //
  // One difficulty is transforming data from the GSL to the LAPACKE
  // complex representations. Luckily, because complex structs are
  // prety similar in general, LAPACKE gives us a way to define their
  // complex type as the gsl_complex type.
  //
  // The second difficulty is that our data are stored in terms of GSL
  // matrices and vectors. These store their data in a row-major order
  // array. LAPACKE routines require us to pass a pointer to this data
  // directly. []->data is that pointer

  char uOption = 'N'; // We don't care about the left eigenvalues.
  char vtOption = 'A'; // We want all of the right eigenvalues.

  gsl_vector* S = gsl_vector_calloc(numBonds);
  gsl_matrix_complex* U 
    = gsl_matrix_complex_calloc(numBonds,numBonds);
  gsl_matrix_complex* VT
    = gsl_matrix_complex_calloc(numBonds,numBonds);
  
  gsl_complex* work 
    = (gsl_complex*) malloc(sizeof(gsl_complex));
  int workLength = -1; // This is a flag to optimize the size
  double* rwork = (double*) malloc(5 * numBonds * sizeof(double));
  
  // Because workLength is -1, the function computes an optimal
  // workspace size, then sets the first element of the work vector to
  // that.
  LAPACKE_zgesvd_work(LAPACK_ROW_MAJOR,
                      uOption, vtOption, 
                      numBonds, numBonds,
                      (gsl_complex*) D->data, numBonds, 
                      S->data, 
                      (gsl_complex*) U->data, numBonds,
                      (gsl_complex*) VT->data, numBonds, 
                      work, workLength, 
                      rwork);
  
  // The first double in the work array is the real part of the first
  // complex number. Since redefining this as an int always rounds
  // down, we need to add 0.5 to ensure we're getting the number that
  // the query returned.
  workLength = (int) GSL_REAL(work[0]) + 0.5; // Set optimal worklength.
  free(work);
  work = (gsl_complex*) malloc(workLength * sizeof(gsl_complex));

  // Now we take the SVD.
  LAPACKE_zgesvd_work(LAPACK_ROW_MAJOR,
                      uOption, vtOption, 
                      numBonds, numBonds,
                      (gsl_complex*) D->data, numBonds, 
                      S->data, 
                      (gsl_complex*) U->data, numBonds,
                      (gsl_complex*) VT->data, numBonds, 
                      work, workLength, 
                      rwork);

  // The singular value vector is ordered in such a way that the last
  // element is the smallest. Presumably, if the characteristic matrix
  // is singular, this will be close to zero. If it is not, return a
  // warning to std::clog, but continue anyway.
  int lastRow = numBonds-1;
  if (gsl_vector_get(S, lastRow) > SVD_EPSILON)
  {
    std::clog << "Warning: smallest singular value was " 
              << gsl_vector_get(S, lastRow) << " > epsilon "
              << SVD_EPSILON
              << std::endl;
  }

  // The conjugate transpose of the last row of the VT matrix 
  // corresponds to the eigenvector of the matrix that corresponds
  // to the singular value. The function returns this whether the
  // smallest singular value was close to zero or not.
  gsl_vector_complex* Eigenvector = gsl_vector_complex_calloc(numBonds);
  for (unsigned int j=0; j<numBonds; j++)
  {
    gsl_complex VLastj 
      = gsl_matrix_complex_get(VT, lastRow, j);
    gsl_vector_complex_set(Eigenvector, j, 
                           gsl_complex_conjugate(VLastj));
  }

  // Clean Up
  free(work);
  free(rwork);
  gsl_matrix_complex_free(D);
  gsl_matrix_complex_free(U);
  gsl_matrix_complex_free(VT);
  gsl_vector_free(S);

  return Eigenvector;
}



gsl_vector_complex* 
QuantumGraph::eigenvectorCheck(const gsl_complex z, const gsl_vector_complex* V)
const
{
  gsl_vector_complex* resultVec = gsl_vector_complex_calloc(numBonds);

  gsl_matrix_complex* D = gsl_matrix_complex_calloc(numBonds, numBonds);
  gsl_matrix_complex_memcpy(D, SMatrix);

  for (unsigned int j=0; j<numBonds; j++)
  {
    gsl_complex negILj
      = gsl_vector_complex_get(negImaginaryLengthVector, j);
    gsl_complex Tj = gsl_complex_exp(gsl_complex_mul(negILj, z));
    gsl_matrix_complex_set(D, j, j, gsl_complex_negative(Tj));
  }

  gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, D,
                 V, GSL_COMPLEX_ONE, resultVec);

  // Clean up.
  gsl_matrix_complex_free(D);

  return resultVec;
}



// Uses much of the same code as the eigenvector code, in the beginning
// but after the SVD is found, that is used to find the value of the
// characteristic function divide by its derivative.
gsl_complex QuantumGraph::newtonStep(const gsl_complex z) const
{
  gsl_matrix_complex* D = gsl_matrix_complex_calloc(numBonds, numBonds);
  gsl_vector_complex* T = gsl_vector_complex_calloc(numBonds);
  gsl_matrix_complex_memcpy(D, SMatrix);

  for (unsigned int j=0; j<numBonds; j++)
  {
    gsl_complex negILj
      = gsl_vector_complex_get(negImaginaryLengthVector, j);
    gsl_complex Tj = gsl_complex_exp(gsl_complex_mul(negILj, z));
    gsl_matrix_complex_set(D, j, j, gsl_complex_negative(Tj));
    gsl_vector_complex_set(T, j, Tj);
  }
  
  // The GSL library does not offer a complex singular value
  // decomposition routine. So we must rely on LAPACK, which has a
  // standard c routine interface called LAPACKE.
  //
  // One difficulty is transforming data from the GSL to the LAPACKE
  // complex representations. Luckily, because complex structs are
  // prety similar in general, LAPACKE gives us a way to define their
  // complex type as the gsl_complex type.
  //
  // The second difficulty is that our data are stored in terms of GSL
  // matrices and vectors. These store their data in a row-major order
  // array. LAPACKE routines require us to pass a pointer to this data
  // directly. []->data is that pointer

  char uOption = 'A'; // We need U
  char vtOption = 'A'; // We need V

  gsl_vector* S = gsl_vector_calloc(numBonds);
  gsl_matrix_complex* U 
    = gsl_matrix_complex_calloc(numBonds,numBonds);
  gsl_matrix_complex* VT
    = gsl_matrix_complex_calloc(numBonds,numBonds);
  
  gsl_complex* work 
    = (gsl_complex*) malloc(sizeof(gsl_complex));
  int workLength = -1; // This is a flag to optimize the size
  double* rwork = (double*) malloc(5 * numBonds * sizeof(double));
  
  // Because workLength is -1, the function computes an optimal
  // workspace size, then sets the first element of the work vector to
  // that.
  LAPACKE_zgesvd_work(LAPACK_ROW_MAJOR,
                      uOption, vtOption, 
                      numBonds, numBonds,
                      (gsl_complex*) D->data, numBonds, 
                      S->data, 
                      (gsl_complex*) U->data, numBonds,
                      (gsl_complex*) VT->data, numBonds, 
                      work, workLength, 
                      rwork);
  
  // The first double in the work array is the real part of the first
  // complex number. Since redefining this as an int always rounds
  // down, we need to add 0.5 to ensure we're getting the number that
  // the query returned.
  workLength = (int) GSL_REAL(work[0]) + 0.5; // Set optimal worklength.
  free(work);
  work = (gsl_complex*) malloc(workLength * sizeof(gsl_complex));

  // Now we take the SVD.
  LAPACKE_zgesvd_work(LAPACK_ROW_MAJOR,
                      uOption, vtOption, 
                      numBonds, numBonds,
                      (gsl_complex*) D->data, numBonds, 
                      S->data, 
                      (gsl_complex*) U->data, numBonds,
                      (gsl_complex*) VT->data, numBonds, 
                      work, workLength, 
                      rwork);

  free(work);
  free(rwork);
  // End SVD work

  // T is now equal to LT.
  gsl_vector_complex_mul(T, LengthVector);

  // Now we find the trace of D^{-1}D' using the SVD, the length vector,
  // and the T vector.
  gsl_complex trace = GSL_COMPLEX_ZERO;
  gsl_complex LTn, Unm, VTmn, mnProduct; 
  double Sm;
  for (unsigned int m = 0; m<numBonds; m++)
  {
    Sm = gsl_vector_get(S, m);
    for (unsigned int n = 0; n<numBonds; n++)
    {
      LTn = gsl_vector_complex_get(T, n);
      Unm = gsl_matrix_complex_get(U, n, m);
      VTmn = gsl_matrix_complex_get(VT, m, n);

      mnProduct = gsl_complex_conjugate(gsl_complex_mul(Unm,VTmn));
      mnProduct = gsl_complex_div_real(mnProduct, Sm);
      mnProduct = gsl_complex_mul(mnProduct, LTn);
      trace = gsl_complex_add(trace, mnProduct);
    }
  }

  trace = gsl_complex_mul_imag(trace, 1.0);

  // Finally, clean up the memory we allocated.
  gsl_matrix_complex_free(D);
  gsl_vector_complex_free(T);
  gsl_matrix_complex_free(U);
  gsl_matrix_complex_free(VT);
  gsl_vector_free(S);

  return gsl_complex_sub(z, gsl_complex_inverse(trace));
}


////////////////////////////////////////////////////////////////////////
// File Output Method


std::ostream & operator<<(std::ostream& os, const QuantumGraph& QG)
{
  // Lose no precision due to storing double as text.
  //os.precision( std::numeric_limits<double>::max_digits10 );
  os.precision(3 );
  os << std::scientific;
  os << std::showpos;

  // Some header information so we know what we're looking at.
  os << "# BEGIN"                      << std::endl;
  os << "#  Quantum Graph Information" << std::endl;
  os << "#  "                          << std::endl;
  os << "#    Number of Bonds:"        << std::endl;
  os << "#      " << QG.numBonds       << std::endl;

  // Write out real part of S-Matrix
  os << "#    S-Matrix (Real Part):"   << std::endl;
  for (unsigned int i=0; i<QG.numBonds; i++)
  {
    os << "#    ";
    for (unsigned int j=0; j<QG.numBonds; j++)
    {
      os << "  " << GSL_REAL(gsl_matrix_complex_get(QG.SMatrix, i, j));
    }
    os << std::endl;
  }
  
  // Write out imaginary part of S-Matrix
  os << "#    S-Matrix (Imaginary Part):" << std::endl;
  for (unsigned int i=0; i<QG.numBonds; i++)
  {
    os << "#    ";
    for (unsigned int j=0; j<QG.numBonds; j++)
    {
      os << "  " << GSL_IMAG(gsl_matrix_complex_get(QG.SMatrix, i, j));
    }
    os << std::endl;
  }

  // Write out real part of bond lengths
  os << "#    Lengths (Real Part):" << std::endl;
  os << "#    ";
  for (unsigned int i=0; i<QG.numBonds; i++)
  {
    os << "  " << GSL_REAL(gsl_vector_complex_get(QG.LengthVector, i));
  }
  os << std::endl;

  // Write out imaginary part of lengths (bond loss factors)
  os << "#    Lengths (Imaginary Part):" << std::endl;
  os << "#    ";
  for (unsigned int i=0; i<QG.numBonds; i++)
  {
    os << "  " << GSL_IMAG(gsl_vector_complex_get(QG.LengthVector, i));
  }
  os << std::endl;

  // So we know when we have all the information about the quantum graph
  // when we examine it in the constructor.
  os << "#  "   << std::endl;
  os << "# END";

  // Works like we expect "os << QG << std::endl;" to do.
  return os;
}











