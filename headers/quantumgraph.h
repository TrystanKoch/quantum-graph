/*----------------------------------------------------------------------
   quantumgraph.cpp - Function Declarations for a class representing a 
                      quantum graph that includes loss.
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

//  Provides the classes
//    QuantumGraph
//
//  QuantumGraph
//    represents a Quantum Graph. This implementation assumes nothing
//    about the scattering vertices or the bonds, except the normal
//    constraint that we describe a quantum graph as a directed graph.
//    The scattering matrices are potentially lossy, and are not assumed
//    to be unitary. The transmission line lengths might also contain a
//    complex part, so that loss during transmission can be modeled.

#ifndef QUANTUMGRAPH_H_
#define QUANTUMGRAPH_H_

#include <ostream>
#include <vector>
#include <string>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


// Contains information about the scattering and transmission of waves
// along a quantum graph. The scattering at vertices is described by the
// S matrix and the transmission along bonds is described by the T 
// matrix, a diagonal matrix that depends on the wavenumber k. 
//
// Given an eigenvalue of the quantum graph, k, then ST(k)a(k) = a(k) 
// for some eigenvector a(k). To find the eigenvalues of the quantum
// graph, we must therefore find roots of the characteristic equation
// Det(1-ST(k)). Then we can find a(k) by finding the null space of
// Det(1-ST(k)). This object does both calculations efficiently.
//
// While it is possible to make the characteristic equation entirely
// real or entirely imaginary for lossless graphs, thereby simplifying
// the root-finding, this class does not assume a lossless problem. Due
// to this, we can manipulate the characteristic equation by multiplying
// it by a number that will never be zero: it is possible to greatly
// reduce the number of operations we must complete this way. Below,
// the user can find the characteristic equation we actually find.
//
// This class takes in the basic information about the quantum graph the
// user has defined in some external way, and performs the basic 
// calculations for finding eigenvalues and eigenvectors.
//
// The characteristic equation that this class finds complex zeros of
// is given by (using LaTeX)
//   \det{[ S - T^{-1}(k) ]} = 0
// which will be zero when the characteristic equation described above
// is zero. The class also finds eigenvectors by finding the null space
// of S-T^{-1}(k) using a singular value decomposition.
//
// Examples of use:
//
// -Construction-
//  If you have the scattering matrix S as a gsl_matrix_complex type and
//  the length vector L as a gsl_vector_complex:
//      QuantumGraph QG = QuantumGraph(L, S);
//  If you have, instead used vectors of doubles RealL and ImagL to
//  describe the bond transmission, and vectors of vectors of doubles
//  RealS and ImagS to describe the vertex scattering, you can have
//      QuantumGraph QG = QuantumGraph(RealL, ImagL, RealS, ImagS);
//  To facilitate storage of the graph's information in a file, you can
//  also take in a header string (format below) to make the definition
//      QuantumGraph QG = QuantumGraph(header);
// 
// -Calculation-
//  To get the complex value of the characteristic function for a
//  complex wavenumber z (gsl_complex type), use
//      gsl_complex c = QG.characteristic(z);
//  Once you have found a wavenumber z that is a zero of the 
//  characteristic equation, you can find the associated eigenvector "a"
//  (gsl_vector_complex type) by
//      gsl_vector_complex *a = QG.eigenvector(z);
//  To check that this is really an eigenvector, you can check to see
//  if it is in the null space of S-T^{-1}(z) by using
//      gsl_vector_complex *n = QG.eigenvectorCheck(z, a);
//  which should be close to the zero vector if it is a solution.
//
class QuantumGraph
{
  private:
    // The number of directed bonds in the graph.
    // Used for many iterations over the bonds and memory allocations.
    unsigned int num_bonds_;


    // The scattering matrix describing the vertex scattering of the
    // quantum graph. May be subunitary in general, for lossy vertex
    // conditions. Unitarity implies no loss at the vertices.
    gsl_matrix_complex* scattering_matrix_;


    // A list of electrical lengths for the bonds in the quantum graph.
    // A real length coresponds to lossless bonds, while a complex
    // length coresponds to a uniform loss along the bond.
    gsl_vector_complex* lengths_;


    // A vector corresponding to -i*lengths_.
    // Because these values show up so often and multiplication requires
    // additional work, the class carries this member for efficiency.
    // Keeping this in vector form lets us skip a matrix multiplication.
    gsl_vector_complex* negative_imaginary_lengths_;


    // The sum of the lengths in lengths_, multiplied by i. Again,
    // This data member exists because the quantity shows up so often
    // that its precalculation saves a noticeable amount of time even
    // though it is not the bottleneck.
    gsl_complex positive_imaginary_lengths_sum_;
    

    // Checks that a matrix has the dimensions of the number of bonds.
    bool CheckMatrixSize(const std::vector<std::vector<double>>) const;


    // Creates negative_imaginary_lengths_ and positive_imaginary_lengths_sum_
    // from the input data. Should be called in every constructor.
    void MakeInternals();


    // These two make the two secondary vectors that cut down on the 
    // number of operations. 
    void MakeNegImaginaryLengthVector();
    void MakePosImaginaryLengthVectorSum();


    // Allocates memory for the GSL structs the class stores.
    void AllocateGraphMemory(const unsigned int);


    // Frees the memory for all of the data members. Invoked by the 
    // destructor to prevent leaking allocated memory.
    void FreeGraphMemory();



  protected:
    // This method is essentially the same as the gsl constructor,
    // Except that it deletes the matrix as it exists and resets it to
    // one described by the new length vector and scattering matrix.
    //
    // This method is intended to let a derived class set the 
    // mathematical portion of the model to match its implementation.
    void Update(const gsl_vector_complex*, const gsl_matrix_complex*);



  public:
    // Default Constructor. Sets the S Matrix and Length Vector to a
    // single entry which is zero.
    QuantumGraph();

    // Basic Constructor. Takes in a list of lengths and a scattering
    // matrix that have already been put into gsl_complex form.
    QuantumGraph(const gsl_vector_complex*, const gsl_matrix_complex*);


    // Vector Constructor. Because can be easiest to store information 
    // in arrays of doubles, this is a way to initialize the object from
    // two vectors (the real and imaginary parts of the bond lengths)
    // and two vectors of vectors (the real and imaginary parts of the
    // Graph scattering matrix).
    QuantumGraph(const std::vector<double>, 
                 const std::vector<double>,
                 const std::vector<std::vector<double>>,
                 const std::vector<std::vector<double>>);


    // String constructor. Because we pass information about the graph
    // to the data files as we're writing them, so that we can keep
    // the problem we're solving consistant across files, we need a way
    // to take the text from a data file's header and initialize our
    // Quantum Graph object based on that data.
    QuantumGraph(const std::vector<std::string>);


    // Standard copy constructor.
    QuantumGraph(const QuantumGraph&);


    // Because there's allocated memory in the object, we have to
    // explicitly free it when we delete the graph. So we have an
    // explicitly defined destructor method.
    ~QuantumGraph();


    // Returns the characteristic function for the quantum graph at a
    // given complex value. When this function returns zero, the 
    // argument corresponds to an eigenvalue of the quantum graph.
    //
    // Input:  gsl_complex -> A complex wavenumber.
    // Output: gsl_complex -> The characteristic equation's value at
    //                        that wavenumber.
    //
    // This finds Det[T^{-1}(z)-S].
    gsl_complex Characteristic(const gsl_complex) const;


    // Finds the eigenvector corresponding to the complex wavenumber z.
    // The vector returned is normalized to a magnitude of 1, which
    // does not correspond to a physical eigenfunction's true elements.
    // However, each eigenvector itself will correspond to an 
    // eigenfunction of the graph.
    //
    // Input:  gsl_complex -> A complex wavenumber (an eigenvalue).
    // Output: gsl_vector_complex* -> The vector corresponding to the
    //                        smallest singular value of the 
    //                        characteristic matrix. (If the input was
    //                        an eigenvalue, this is an eigenvector of
    //                        the quantum graph.)
    //
    // While this is a step in finding normalized eigenfunctions, the
    // eventual normalization will rely on a physical analog to a
    // particular problem.
    gsl_vector_complex* Eigenvector(const gsl_complex) const;


    // Checks the eigenvector by multiplying it by the T^-1 - S 
    // characteristic matrix. If the resulting vector is very close to
    // the zero vector, then we can trust that the eigenvector is really
    // an eigenvector of the quantum graph we're studying. Note that
    // because the characteristic equation changes depending on the
    // eigenvalue, we need it as an input here.
    //
    // Input:  gsl_complex -> A complex wavenumber (an eigenvalue).
    // Input:  gsl_vector_complex* -> A potential eigenvector.
    // Output: double -> The complex vector 2-norm of the characteristic 
    //                   matrix times the vector you are testing. For a
    //                   true eigenvalue, eigenvector pair, this will be
    //                   very close to zero.
    //
    // This is less for computing things of interest, and more of a
    // check that all the above steps worked out when finding the
    // eigenvector.
    double EigenvectorCheck(const gsl_complex, 
                            const gsl_vector_complex*) const;


    // Given an initial guess z, computes z-f(z)/f'(z), to use in
    // Newton's root finding method. This method uses the singluar
    // value decomposition to find an inverse, then uses that
    // and other quantities already found above to calculate the
    // determinant
    //
    // Input:  gsl_complex -> A complex value to start the current 
    //                        Newton's method iteration.
    // Output: gsl_complex -> The starting point for the next Newton's
    //                        method iteration.
    //
    // This only computes one step in the iterative root-finding method.
    // It must be called by a loop that performs the iteration. This
    // method is a part of the class because it is easier to perform the
    // singular value decomposition with access to the private data
    // members.
    gsl_complex NewtonStep(const gsl_complex) const;
    
    
    // Given a Quantum Graph object with suitable definitions for the
    // bond lengths, this function normalizes that graph, keeping the
    // bonds proportional while setting the total length so that the
    // inter-eigenvalue spacing is normalized. It sums up the bond
    // lengths then multiplies each length by 2*pi and divides by the
    // length sum. This essentially scales the bond lengths such that
    // the total length of the undirected bonds is pi. From weyl's
    // formula, we can see that the average eigenvalue spacing is one.
    //
    // This method changes the object! It uses the update method above.
    // Nothing is returned, but the original quantum graph is lost and
    // replaced.
    void Normalize();


    // Prints the information about a quantum graph into an output 
    // stream for clearer file I/O. The real and imaginary parts of the
    // scattering matrix and the length vector are output using "#" as
    // the first character in the line. This is intended to serve as a
    // comment that gnuplot and other programs can ignore, while still
    // containing information about what problem the computed data comes
    // from. 
    //
    // Input:  QuantumGraph& -> The quantum graph object.
    // Output: std::ostream& -> An output stream full of information
    //                          about the quantum Graph.
    //
    // The text output from this function is what the string-based
    // constructor above will use to create a copy of the orginal graph.
    // this way, we can pass the same problem between different programs
    // without having to manually keep track of what files match which
    // problem.
    friend std::ostream& operator<<(std::ostream&, 
                                    const QuantumGraph&);
};

#endif
