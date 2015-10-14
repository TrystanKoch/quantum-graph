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
  num_bonds_ = 1;

  AllocateGraphMemory(num_bonds_);

  MakeInternals();
}

// Constructor from GSL types
QuantumGraph::QuantumGraph(const gsl_vector_complex* lengths,
                           const gsl_matrix_complex* scattering_matrix)
{
  num_bonds_ = scattering_matrix->size1;

  if (lengths->size != num_bonds_)
  {
    std::clog << "Length Vector and S Matrix must be same dimension." 
              << std::endl;
  }
  else if (scattering_matrix->size2 != num_bonds_)
  {
    std::clog << "S Matrix must be square." 
              << std::endl;
  }

  AllocateGraphMemory(num_bonds_);

  gsl_matrix_complex_memcpy(scattering_matrix_, scattering_matrix);
  gsl_vector_complex_memcpy(lengths_, lengths);

  MakeInternals();
}




// Constructor from vectors of doubles representing real and imaginary
// parts of lengths_ and scattering_matrix_.
QuantumGraph::QuantumGraph(
    const std::vector<double> lengths_real, 
    const std::vector<double> lengths_imag,
    const std::vector<std::vector<double>> scattering_matrix_real,
    const std::vector<std::vector<double>> scattering_matrix_imag
)
{
  // Make sure we're not giving the object vectors and matrices of
  // different sizes.
  num_bonds_ = lengths_real.size();
  if (num_bonds_ != lengths_imag.size())
  {
    std::clog << "Imaginary part of Length Vector has incorrect size." 
              << std::endl;
  }
  if (!CheckMatrixSize(scattering_matrix_real))
  {
    std::clog << "Imaginary part of S Matrix has incorrect size." 
              << std::endl;
  }
  if (!CheckMatrixSize(scattering_matrix_imag))
  {
    std::clog << "Real part of S Matrix has incorrect size." 
              << std::endl;
  }

  // Okay, allocate the memory for data members.
  AllocateGraphMemory(num_bonds_);
  
  // Make scattering_matrix_
  gsl_complex Sij;
  for (unsigned int i=0; i<num_bonds_; i++)
  {
    for (unsigned int j=0; j<num_bonds_; j++)
    {
      Sij = gsl_complex_rect(scattering_matrix_real[i][j],
                             scattering_matrix_imag[i][j]);
      gsl_matrix_complex_set(scattering_matrix_, i, j, Sij);
    }
  }

  // Make lengths_
  for (unsigned int j=0; j<num_bonds_; j++)
  {
    gsl_complex Lj 
        = gsl_complex_rect(lengths_real[j], lengths_imag[j]);
    gsl_vector_complex_set(lengths_, j, Lj);
  }

  MakeInternals();
}


// Constructor that examines a file header.
//  Each program's output has a header that contains the details of the
//  quantum graph. This constructor is for reading a string containing
//  that header and initializing a Quantum Graph based on it.
//  See the header file for information about the formatting.
QuantumGraph::QuantumGraph(std::vector<std::string> header)
{
  std::vector<std::vector<double>> scattering_matrix_real;
  std::vector<std::vector<double>> scattering_matrix_imag;
  std::vector<double> lengths_imag;
  std::vector<double> lengths_real;

  num_bonds_ = 0;

  for (unsigned int i=0; i<header.size(); i++)
  {
    std::string line = header[i];
    if (line == "#    Number of Bonds:")
    {
      line = header[i+1];
      line.erase(0,1);
      num_bonds_ = std::stoi(line);
      lengths_imag.resize(num_bonds_);
      lengths_real.resize(num_bonds_);
      scattering_matrix_imag.resize(num_bonds_);
      scattering_matrix_real.resize(num_bonds_);
      for (unsigned int m=0; m<num_bonds_; m++)
      {
        scattering_matrix_imag[m].resize(num_bonds_);
        scattering_matrix_real[m].resize(num_bonds_);
      }
    }
    else if (line == "#    S-Matrix (Real Part):")
    {
      for (unsigned int m=0; m<num_bonds_; m++)
      {
        line = header[i+1+m];
        line.erase(0,1);
        std::stringstream scattering_matrix_real_stream(line);
        for (unsigned n=0; n<num_bonds_; n++)
        {
          scattering_matrix_real_stream >> scattering_matrix_real[m][n];
        }
      }
    }
    else if (line == "#    S-Matrix (Imaginary Part):")
    {
      for (unsigned int m=0; m<num_bonds_; m++)
      {
        line = header[i+m+1];
        line.erase(0,1);
        std::stringstream scattering_matrix_imag_stream(line);

        for (unsigned n=0; n<num_bonds_; n++)
        {
          scattering_matrix_imag_stream >> scattering_matrix_imag[m][n];
        }
      }
    }
    else if (line == "#    Lengths (Real Part):")
    {
      line = header[i+1];
      line.erase(0,1);
      std::stringstream lengths_real_stream(line);
      for (unsigned int m=0; m<num_bonds_; m++)
      {
        lengths_real_stream >> lengths_real[m];
      }
    }
    else if (line == "#    Lengths (Imaginary Part):")
    {
      line = header[i+1];
      line.erase(0,1);
      std::stringstream lengths_imag_stream(line);
      for (unsigned int m=0; m<num_bonds_; m++)
      {
        lengths_imag_stream >> lengths_imag[m];
      }
    }
    else
    {
      continue;
    }
  }

  if (num_bonds_!=lengths_real.size())
  {
    std::clog << "Real part of Length Vector has incorrect size." 
              << std::endl;
  }
  if (num_bonds_!=lengths_imag.size())
  {
    std::clog << "Imaginary part of Length Vector has incorrect size." 
              << std::endl;
  }
  if (!CheckMatrixSize(scattering_matrix_real))
  {
    std::clog << "Imaginary part of S Matrix has incorrect size." 
              << std::endl;
  }
  if (!CheckMatrixSize(scattering_matrix_imag))
  {
    std::clog << "Real part of S Matrix has incorrect size." 
              << std::endl;
  }
  

  AllocateGraphMemory(num_bonds_);

  for (unsigned int i=0; i<num_bonds_; i++)
  {
    for (unsigned int j=0; j<num_bonds_; j++)
    {
      gsl_complex Sij = gsl_complex_rect(scattering_matrix_real[i][j],
                                         scattering_matrix_imag[i][j]);
      gsl_matrix_complex_set(scattering_matrix_, i, j, Sij);
    }
  }

  for (unsigned int i=0; i<num_bonds_; i++)
  {
    gsl_complex Li = gsl_complex_rect(lengths_real[i],
                                      lengths_imag[i]);
    gsl_vector_complex_set(lengths_, i, Li);
  }

  MakeInternals();
}




// Copy constructor.
QuantumGraph::QuantumGraph(const QuantumGraph& QG)
{
  num_bonds_ = QG.num_bonds_;
  AllocateGraphMemory(num_bonds_);

  gsl_matrix_complex_memcpy(scattering_matrix_, QG.scattering_matrix_);
  gsl_vector_complex_memcpy(lengths_, QG.lengths_);

  MakeInternals();
}




// Need to get rid of allocated data when we destroy the object, because
// GSL uses pure C and needs to allocate memory.
QuantumGraph::~QuantumGraph()
{
  FreeGraphMemory();
}




////////////////////////////////////////////////////////////////////////
// Internal Data initialization Methods


void QuantumGraph::MakeInternals()
{
  MakeNegImaginaryLengthVector();
  MakePosImaginaryLengthVectorSum();
}




// So we don't have to do this on every computation.
void QuantumGraph::MakeNegImaginaryLengthVector()
{
  gsl_vector_complex_memcpy(negative_imaginary_lengths_, lengths_);
  gsl_vector_complex_scale(negative_imaginary_lengths_, 
                           gsl_complex_rect(0, -1) );
}




// So we don't have to do this on every computation.
void QuantumGraph::MakePosImaginaryLengthVectorSum()
{
  gsl_complex length_sum = GSL_COMPLEX_ZERO;
  for (unsigned int i=0; i<num_bonds_; i++)
  {
    gsl_complex Li = gsl_vector_complex_get(lengths_, i);
    length_sum = gsl_complex_add(length_sum, Li);
  }
  positive_imaginary_lengths_sum_ 
    = gsl_complex_mul(length_sum, gsl_complex_rect(0, 1));
}


// Check the size of a vector of vectors during initialization.
bool QuantumGraph::CheckMatrixSize(
    const std::vector<std::vector<double>> matrix
) const
{
  if (matrix.size() != num_bonds_)
  {
    return 0;
  }
  for (unsigned int i=0; i<num_bonds_; i++)
  {
    if (matrix[i].size() != num_bonds_)
    {
      return 0;
    }
  }
  return 1;
}




////////////////////////////////////////////////////////////////////////
// Memory Handling Methods


void QuantumGraph::AllocateGraphMemory(const unsigned int N)
{
  scattering_matrix_ = gsl_matrix_complex_calloc(N,N);
  lengths_ = gsl_vector_complex_calloc(N);
  negative_imaginary_lengths_ = gsl_vector_complex_calloc(N);
}




void QuantumGraph::FreeGraphMemory()
{
  gsl_matrix_complex_free(scattering_matrix_);
  gsl_vector_complex_free(lengths_);
  gsl_vector_complex_free(negative_imaginary_lengths_);
}

////////////////////////////////////////////////////////////////////////
/// Methods for Derived Classes only

// Protected method for resetting the representation of the graph based
// on the implementation of a derived class.
void QuantumGraph::Update(const gsl_vector_complex* lengths,
                          const gsl_matrix_complex* scattering_matrix)
{
  FreeGraphMemory();

  num_bonds_ = scattering_matrix->size1;

  if (lengths->size != num_bonds_)
  {
    std::clog << "Length Vector and S Matrix must be same dimension." 
              << std::endl;
  }
  else if (scattering_matrix->size2 != num_bonds_)
  {
    std::clog << "S Matrix must be square." 
              << std::endl;
  }

  AllocateGraphMemory(num_bonds_);

  gsl_matrix_complex_memcpy(scattering_matrix_, scattering_matrix);
  gsl_vector_complex_memcpy(lengths_, lengths);

  MakeInternals();
}


////////////////////////////////////////////////////////////////////////
// Computation Methods


// While this function will have zeros, it will not be strictly real
// along the real axis. Only a complex zero-finding algorithm will work.
gsl_complex QuantumGraph::Characteristic(const gsl_complex z) const
{
  gsl_matrix_complex* characteristic_matrix 
      = gsl_matrix_complex_calloc(num_bonds_, num_bonds_);
  gsl_matrix_complex_memcpy(characteristic_matrix, scattering_matrix_);

  for (unsigned int j=0; j<num_bonds_; j++)
  {
    gsl_complex neg_imag_Lj
      = gsl_vector_complex_get(negative_imaginary_lengths_, j);
    gsl_complex Tj = gsl_complex_exp(gsl_complex_mul(neg_imag_Lj, z));
    gsl_matrix_complex_set(characteristic_matrix, j, j, 
                           gsl_complex_negative(Tj));
  }

  int signum;
  gsl_permutation* permutation = gsl_permutation_alloc(num_bonds_);

  gsl_linalg_complex_LU_decomp(characteristic_matrix, 
                               permutation, &signum);
  gsl_complex determinant 
      = gsl_linalg_complex_LU_det(characteristic_matrix, signum);

  gsl_permutation_free(permutation);
  gsl_matrix_complex_free(characteristic_matrix);

  return determinant;
}





gsl_vector_complex* QuantumGraph::Eigenvector(const gsl_complex z) const
{
  // For simplicity, D will represent the characteristic matrix at a
  // given input z. 
  gsl_matrix_complex* D 
    = gsl_matrix_complex_calloc(num_bonds_, num_bonds_);
  gsl_matrix_complex_memcpy(D, scattering_matrix_);

  for (unsigned int j=0; j<num_bonds_; j++)
  {
    gsl_complex neg_imag_Lj
      = gsl_vector_complex_get(negative_imaginary_lengths_, j);
    gsl_complex Tj = gsl_complex_exp(gsl_complex_mul(neg_imag_Lj, z));
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

  char u_option = 'N'; // We don't care about the left eigenvalues.
  char vt_option = 'A'; // We want all of the right eigenvalues.

  // The variables are based off of
  //   D = USV^H
  // where H is the hermitian conjugate, and the other three are the
  // result of singular value decomposition.
  gsl_vector* S = gsl_vector_calloc(num_bonds_);
  gsl_matrix_complex* U 
    = gsl_matrix_complex_calloc(num_bonds_,num_bonds_);
  gsl_matrix_complex* VT
    = gsl_matrix_complex_calloc(num_bonds_,num_bonds_);
  
  gsl_complex* work 
    = (gsl_complex*) malloc(sizeof(gsl_complex));
  int workLength = -1; // This is a flag to optimize the size
  double* rwork = (double*) malloc(5 * num_bonds_ * sizeof(double));
  
  // Because workLength is -1, the function computes an optimal
  // workspace size, then sets the first element of the work vector to
  // that.
  LAPACKE_zgesvd_work(LAPACK_ROW_MAJOR,
                      u_option, vt_option, 
                      num_bonds_, num_bonds_,
                      (gsl_complex*) D->data, num_bonds_, 
                      S->data, 
                      (gsl_complex*) U->data, num_bonds_,
                      (gsl_complex*) VT->data, num_bonds_, 
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
                      u_option, vt_option, 
                      num_bonds_, num_bonds_,
                      (gsl_complex*) D->data, num_bonds_, 
                      S->data, 
                      (gsl_complex*) U->data, num_bonds_,
                      (gsl_complex*) VT->data, num_bonds_, 
                      work, workLength, 
                      rwork);

  // The singular value vector is ordered in such a way that the last
  // element is the smallest. Presumably, if the characteristic matrix
  // is singular, this will be close to zero. If it is not, return a
  // warning to std::clog, but continue anyway.
  int lastRow = num_bonds_ - 1;
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
  gsl_vector_complex* eigenvector = gsl_vector_complex_calloc(num_bonds_);
  for (unsigned int j=0; j<num_bonds_; j++)
  {
    gsl_complex VLastj 
      = gsl_matrix_complex_get(VT, lastRow, j);
    gsl_vector_complex_set(eigenvector, j, 
                           gsl_complex_conjugate(VLastj));
  }

  // Clean Up
  free(work);
  free(rwork);
  gsl_matrix_complex_free(D);
  gsl_matrix_complex_free(U);
  gsl_matrix_complex_free(VT);
  gsl_vector_free(S);

  return eigenvector;
}



double
QuantumGraph::EigenvectorCheck(const gsl_complex z, 
                               const gsl_vector_complex* v)
const
{
  gsl_vector_complex* result = gsl_vector_complex_calloc(num_bonds_);

  gsl_matrix_complex* characteristic_matrix
      = gsl_matrix_complex_calloc(num_bonds_, num_bonds_);
  gsl_matrix_complex_memcpy(characteristic_matrix, scattering_matrix_);

  for (unsigned int j=0; j<num_bonds_; j++)
  {
    gsl_complex neg_imag_Lj
      = gsl_vector_complex_get(negative_imaginary_lengths_, j);
    gsl_complex Tj = gsl_complex_exp(gsl_complex_mul(neg_imag_Lj, z));
    gsl_matrix_complex_set(characteristic_matrix, j, j, 
                           gsl_complex_negative(Tj));
  }

  gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, characteristic_matrix,
                 v, GSL_COMPLEX_ONE, result);

  // Clean up.
  gsl_matrix_complex_free(characteristic_matrix);

  return gsl_blas_dznrm2(result);
}



// Uses much of the same code as the eigenvector code, in the beginning
// but after the SVD is found, that is used to find the value of the
// characteristic function divide by its derivative.
gsl_complex QuantumGraph::NewtonStep(const gsl_complex z) const
{
  // Once again, D will stand for our characteristic matrix
  // T will stand for our transmission matrix, which we need seperately 
  // this time (before, we'd only need it as part of the characteristic
  // matrix, but now it comes in at a different place in our derivative.
  gsl_matrix_complex* D = gsl_matrix_complex_calloc(num_bonds_, num_bonds_);
  gsl_vector_complex* T = gsl_vector_complex_calloc(num_bonds_);
  gsl_matrix_complex_memcpy(D, scattering_matrix_);

  for (unsigned int j=0; j<num_bonds_; j++)
  {
    gsl_complex negILj
      = gsl_vector_complex_get(negative_imaginary_lengths_, j);
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
  char u_option = 'A'; // We need U
  char vt_option = 'A'; // We need V

  // The variables are based off of
  //   D = USV^H
  // where H is the hermitian conjugate, and the other three are the
  // result of singular value decomposition.
  gsl_vector* S = gsl_vector_calloc(num_bonds_);
  gsl_matrix_complex* U 
    = gsl_matrix_complex_calloc(num_bonds_,num_bonds_);
  gsl_matrix_complex* VT
    = gsl_matrix_complex_calloc(num_bonds_,num_bonds_);
  
  gsl_complex* work 
    = (gsl_complex*) malloc(sizeof(gsl_complex));
  int workLength = -1; // This is a flag to optimize the size
  double* rwork = (double*) malloc(5 * num_bonds_ * sizeof(double));
  
  // Because workLength is -1, the function computes an optimal
  // workspace size, then sets the first element of the work vector to
  // that.
  LAPACKE_zgesvd_work(LAPACK_ROW_MAJOR,
                      u_option, vt_option, 
                      num_bonds_, num_bonds_,
                      (gsl_complex*) D->data, num_bonds_, 
                      S->data, 
                      (gsl_complex*) U->data, num_bonds_,
                      (gsl_complex*) VT->data, num_bonds_, 
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
                      u_option, vt_option, 
                      num_bonds_, num_bonds_,
                      (gsl_complex*) D->data, num_bonds_, 
                      S->data, 
                      (gsl_complex*) U->data, num_bonds_,
                      (gsl_complex*) VT->data, num_bonds_, 
                      work, workLength, 
                      rwork);

  free(work);
  free(rwork);
  // End SVD work

  // T is now equal to L*T.
  gsl_vector_complex_mul(T, lengths_);

  // Now we find the trace of D^{-1}D' using the SVD, the length vector,
  // and the T vector.
  gsl_complex trace = GSL_COMPLEX_ZERO;
  gsl_complex LTn, Unm, VTmn, mnProduct; 
  double Sm;
  for (unsigned int m = 0; m<num_bonds_; m++)
  {
    Sm = gsl_vector_get(S, m);
    for (unsigned int n = 0; n<num_bonds_; n++)
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

  // Finally, return the next iteration for the Newton Root finder.
  // This is our next guess for where the root is.
  return gsl_complex_sub(z, gsl_complex_inverse(trace));
}


////////////////////////////////////////////////////////////////////////
// File Output Method


std::ostream & operator<<(std::ostream& os, const QuantumGraph& QG)
{
  // Lose no precision due to storing double as text.
  os.precision( std::numeric_limits<double>::max_digits10 );
  //os.precision(3); //debuging
  os << std::scientific;
  os << std::showpos;

  // Some header information so we know what we're looking at.
  os << "# BEGIN"                      << std::endl;
  os << "#  Quantum Graph Information" << std::endl;
  os << "#  "                          << std::endl;
  os << "#    Number of Bonds:"        << std::endl;
  os << "#      " << QG.num_bonds_       << std::endl;

  // Write out real part of S-Matrix
  os << "#    S-Matrix (Real Part):"   << std::endl;
  for (unsigned int i=0; i<QG.num_bonds_; i++)
  {
    os << "#    ";
    for (unsigned int j=0; j<QG.num_bonds_; j++)
    {
      os << "  " 
         << GSL_REAL(gsl_matrix_complex_get(QG.scattering_matrix_,i,j));
    }
    os << std::endl;
  }
  
  // Write out imaginary part of S-Matrix
  os << "#    S-Matrix (Imaginary Part):" << std::endl;
  for (unsigned int i=0; i<QG.num_bonds_; i++)
  {
    os << "#    ";
    for (unsigned int j=0; j<QG.num_bonds_; j++)
    {
      os << "  " 
         << GSL_IMAG(gsl_matrix_complex_get(QG.scattering_matrix_,i,j));
    }
    os << std::endl;
  }

  // Write out real part of bond lengths
  os << "#    Lengths (Real Part):" << std::endl;
  os << "#    ";
  for (unsigned int i=0; i<QG.num_bonds_; i++)
  {
    os << "  " << GSL_REAL(gsl_vector_complex_get(QG.lengths_, i));
  }
  os << std::endl;

  // Write out imaginary part of lengths (bond loss factors)
  os << "#    Lengths (Imaginary Part):" << std::endl;
  os << "#    ";
  for (unsigned int i=0; i<QG.num_bonds_; i++)
  {
    os << "  " << GSL_IMAG(gsl_vector_complex_get(QG.lengths_, i));
  }
  os << std::endl;

  // So we know when we have all the information about the quantum graph
  // when we examine it in the constructor.
  os << "#  "   << std::endl;
  os << "# END";

  // Works like we expect "os << QG << std::endl;" to do.
  return os;
}











