/*----------------------------------------------------------------------
   Perturb_quotient_graph.cpp 
     - Choose a doublet in the GSE symmetric graph, perturb a bond to
       see how the symmetry breaks.
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



#include "quantumgraphobject.h"
#include "quantumgraph.h"
#include "quantumgraphrootfinding.h"

#include <cmath>
#include <limits>
#include <string>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_qrng.h>

#include <iostream>
#include <sstream>
#include <streambuf>


// These constants define the boundary of the box in the complex plane.
// X denotes the real line, Y denotes the imaginary line.
//
// Right now this is the only way to specify a particular box. This may 
// change in the future if I figure out an elegent way to pass these four
// numbers in a way that is understandable.
const double X_MIN = +10;
const double X_MAX = +14;
const double Y_MIN = -0.5;
const double Y_MAX = +0.51;


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
const unsigned int NUM_HALVINGS_X = 10;
const unsigned int NUM_HALVINGS_Y = 10;


// Function that returns a graph that is ever so slightly changed in one
// bond of the graph. In this program, we perturb the GSE graph.
QuantumGraph PerturbGraph(const double);

// Computes the splitting for a small perturbation.
double Splitting(const gsl_complex, const gsl_complex);

// This will govern how small our perturbation steps are.
const double STEP_SIZE = 1E-4;
const int NUM_STEPS = 5;

const unsigned int REAL_HALVES = 10;
const unsigned int IMAG_HALVES = 10;

int main()
{
  
  // This will be our initial unperturbed graph
  QuantumGraph QGO;
  QGO = PerturbGraph(0);
  


  // Since I have foolishly written the root bounding and root finding
  // algorithms to output to cout (I wanted to make things modular), and
  // now I want them contained within one program, I need to steal cout
  // and redirect that to something I can read here.
  std::streambuf* old_cout_stream_buffer = std::cout.rdbuf();
  std::stringstream str_cout;
  std::cout.rdbuf(str_cout.rdbuf());
  
  std::clog << "Finding Root Bounding Boxes...  " << std::flush;
  // Now we bound the roots for the original graph.
  BoundingBox initialBoundary;
  initialBoundary 
            = initializeBoundingBox(X_MIN, X_MAX, Y_MIN, Y_MAX, 
                                    NUM_HALVINGS_X, NUM_HALVINGS_Y, QGO);
  std::clog << "Done" << std::endl;
                                    
  std::clog << "Winding Number is " << windingNumber(initialBoundary) << std::endl;
  
  unsigned long long int total;
  total = recursiveRootBounding(initialBoundary, QGO, 2);
  std::clog<< "There are "<< total << " roots."<< std::endl;
  
  
  std::clog << "Creating Bounding Boxes...  " << std::endl;
  // Using our string stream, we create our vector of bounding boxes
  std::stringstream ss;
  std::vector<BoundingBox> bounding_boxes;
  std::string line;
  while (std::getline(str_cout, line))
  {
    BoundingBox BB;
    
    ss.str(line);
    ss >> BB.realMin;
    ss >> BB.realMax;
    ss >> BB.imagMin;
    ss >> BB.imagMax;
    ss.clear();
    
    bounding_boxes.push_back(BB);
    
    //debug
    std::clog << BB <<std::endl;
  }
  std::clog << "There are "<< bounding_boxes.size() 
            << " root doublets."<< std::endl;
 // str_cout.clear();
            
            


  
  
 
  
  // Now we find those double roots and place them in an array
  std::vector<gsl_complex> roots;
  gsl_complex root;
  gsl_qrng* qrng = gsl_qrng_alloc(gsl_qrng_sobol, 2);
  for (unsigned int j=0; j<bounding_boxes.size(); j++)
  {
    std::clog << "\rFinding Roots...  " << j << std::flush;
    BoundingBox BB = bounding_boxes[j];
    root = newtonRootFinderQRNG(BB, qrng, QGO);
    roots.push_back(root);
  }
  //std::cout << std::endl;
  std::clog << "\rFinding Roots...  Done" << std::endl;
  


  
  // We find the second derivative of the secular function at each
  // root.
 // std::vector<gsl_complex> second_derivatives;
 // gsl_complex second_derivative;
 // for (unsigned int i=0; i<roots.size(); i++)
 // {
 //   std::clog << "\rFinding Second Derivatives... " << i << std::flush;
 //   second_derivative = QGO.CharacteristicDerivative(roots[i], 2);
 //   second_derivatives.push_back(second_derivative);
 // }
 // std::clog << "\rFinding Second Derivatives...  Done" << std::endl;
  
  std::cout << "1"<<std::endl;

  // Now we perturb the graph, and for each root, we evaluate the
  // perturbed secular equation at the old root, and use that number
  // to determine the spacing between the two now non-degenerate roots.
  QuantumGraph QG;
 // gsl_complex f, sd;
  std::vector<double> splittings;
  BoundingBox BB;
  for (int i=1; i<NUM_STEPS; i++)
  {
    //std::cout.clear();
   // std::cout.rdbuf(str_cout.rdbuf());
    splittings.erase(splittings.begin(), splittings.end());
    //std::clog << "Perturbing Quantum Graph...  " << i << std::endl;
    std::cout << "2"<<std::endl;
    QG = PerturbGraph(STEP_SIZE*i);
    std::clog << "\nPerturbation is " << STEP_SIZE*i<< std::endl;
    BB = bounding_boxes[i-1];
    
    //BB = //FIXME  Boundary needs to be instantiated
    
    for (unsigned int j=0; j<roots.size(); j++)
    {
      std::cout << "3"<<std::endl;
      std::clog << "Root is:  " << GSL_REAL(roots[j]) <<std::endl;
      std::vector<gsl_complex> roots_doublet;
     // f = QG.Characteristic(roots[j]);
     std::cout << "4"<<std::endl;
     // sd = second_derivatives[j];
     // double approx_splitting = 3*Splitting(f, sd);
      
   /*   BB = initializeBoundingBox(GSL_REAL(roots[j]) - approx_splitting,
                                 GSL_REAL(roots[j]) + approx_splitting,
                                 GSL_IMAG(roots[j]) - approx_splitting,
                                 GSL_IMAG(roots[j]) + approx_splitting,
                                 REAL_HALVES, IMAG_HALVES, QG);
                 std::clog<< "Bounding Box is:  "<< BB<< std::endl;*/
                 
     // std::clog << windingNumber(BB) <<std::endl;
         std::clog<< BB<<std::endl; 
         std::cout.rdbuf(old_cout_stream_buffer);
         std::clog<< BB<<std::endl; 
      recursiveRootBounding(BB, QG, 1);
      std::clog<< "got here" <<std::endl;
      std::vector<BoundingBox> temp_bounding_boxes;                  
      while (std::getline(str_cout, line))
      {
        BoundingBox BB_temp;
    
        ss.str(line);
        //std::clog << line << std::endl;
        ss >> BB_temp.realMin;
        ss >> BB_temp.realMax;
        ss >> BB_temp.imagMin;
        ss >> BB_temp.imagMax;
        ss.clear();
    
        temp_bounding_boxes.push_back(BB_temp);
    
        //debug
        std::clog << BB_temp <<std::endl;
      }
      //std::clog << "Finding Roots...  " << j << std::endl;
      for (unsigned int l=0; l<temp_bounding_boxes.size(); l++)
      {
        root = newtonRootFinderQRNG(temp_bounding_boxes[l], qrng, QG);
        roots_doublet.push_back(root);
      }
      std::clog << "("<< GSL_REAL(roots_doublet[0])<< ","<<GSL_REAL(roots_doublet[1])  <<")" <<std::endl;
      //std::clog<< roots_doublet.size()<< std::endl;
      splittings.push_back(gsl_complex_abs(gsl_complex_sub(roots_doublet[0],roots_doublet[1])));
      //std::clog << "\t" << approx_splitting <<std::flush;
      //std::clog << "\t" << splittings[j] <<std::endl;
      //std::clog << "\rPerturbation " << i << ",  Root number " << j << std::flush;
      temp_bounding_boxes.erase(temp_bounding_boxes.begin(), temp_bounding_boxes.end());
      roots_doublet.erase(roots_doublet.begin(), roots_doublet.end());
      std::clog << "Doublet array size is" << roots_doublet.size()<< std::endl;
    }
    
    

    std::cout.rdbuf(old_cout_stream_buffer);
    std::cout << STEP_SIZE*i;
    for (unsigned int j=0; j<splittings.size(); j++)
    {
      std::cout << "\t" << splittings[j];
    }
    std::cout << std::endl;
    std::cout.rdbuf(str_cout.rdbuf());
    
    
  }
    
  std::cout.rdbuf(old_cout_stream_buffer);
  
  std::clog << "\rPerturbing Quantum Graph...  Done" << std::endl;
  std::clog << roots.size() << " Roots Perturbed" << std::endl;


  
  return 0;
}

double Splitting(const gsl_complex f, const gsl_complex sd)
{
  return std::sqrt(8.0*gsl_complex_abs(f)/gsl_complex_abs(sd));
}


QuantumGraph PerturbGraph(const double x)
{
  gsl_complex length1 = gsl_complex_rect(std::sqrt(1), 0);
   
  gsl_complex length2 = gsl_complex_rect(std::sqrt(2), 0);
  
  gsl_complex length3 = gsl_complex_rect(std::sqrt(3), 0);
  
  gsl_complex length4 = gsl_complex_rect(std::sqrt(5), 0);
  
  gsl_complex length5 = gsl_complex_rect(std::sqrt(6), 0);
  
  gsl_complex length6 = gsl_complex_rect(std::sqrt(7), 0);
  
  gsl_complex length7 = gsl_complex_rect(std::sqrt(10), 0);
  
  gsl_complex length8 = gsl_complex_rect(std::sqrt(11), 0);
  
  gsl_complex length9 = gsl_complex_rect(std::sqrt(13), 0);
  
  gsl_complex length10 = gsl_complex_rect(std::sqrt(14), 0);
  
  gsl_complex length11 = gsl_complex_rect(1/std::sqrt(15), 0);
  
  
  /// Define 3-bond Neumann node  
  std::vector<std::vector<double>> N3 = 
  {
    {-1./3.,  2./3.,  2./3.},
    { 2./3., -1./3.,  2./3.},
    { 2./3.,  2./3., -1./3.}
  };
  
  gsl_matrix_complex* neumann3 = gsl_matrix_complex_calloc(3, 3);
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    { 
      gsl_matrix_complex_set(neumann3, i, j, gsl_complex_rect(N3[i][j], 0));
    }
  }
  
  /// Define 4-bond Neumann node  
  std::vector<std::vector<double>> N4 = 
  {
    {-1./2.,  1./2.,  1./2.,  1./2.},
    { 1./2., -1./2.,  1./2.,  1./2.},
    { 1./2.,  1./2., -1./2.,  1./2.},
    { 1./2.,  1./2.,  1./2., -1./2.}
  };

  gsl_matrix_complex* neumann4 = gsl_matrix_complex_calloc(4, 4);
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
    { 
      gsl_matrix_complex_set(neumann4, i, j, gsl_complex_rect(N4[i][j], 0));
    }
  }
  
  /// Define 3-bond circulator node  (left orientation)
  std::vector<std::vector<double>> LC3 = 
  {
    { 0,  1,  0},
    { 0,  0,  1},
    { 1,  0,  0}
  };
  
  gsl_matrix_complex* left_circulator3 = gsl_matrix_complex_calloc(3, 3);
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    { 
      gsl_matrix_complex_set(left_circulator3, i, j, gsl_complex_rect(LC3[i][j], 0));
    }
  }
  
  
  /// Define 3-bond circulator node  (right orientation)
  std::vector<std::vector<double>> RC3 = 
  {
    { 0,  0,  1},
    { 1,  0,  0},
    { 0,  1,  0}
  };
  
  gsl_matrix_complex* right_circulator3 = gsl_matrix_complex_calloc(3, 3);
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    { 
      gsl_matrix_complex_set(right_circulator3, i, j, gsl_complex_rect(RC3[i][j], 0));
    }
  }


  /// Define 1-bond Neumann node (open reflector)
  gsl_matrix_complex* neumann_reflector = gsl_matrix_complex_calloc(1,1);
  gsl_matrix_complex_set(neumann_reflector, 0, 0, GSL_COMPLEX_ONE);

  /// Define 1-bond Dirichlet node (closed reflector)
  gsl_matrix_complex* dirichlet_reflector = gsl_matrix_complex_calloc(1,1);
  gsl_matrix_complex_set(dirichlet_reflector, 0, 0, GSL_COMPLEX_NEGONE);
  
  


 
  gsl_complex length1_perturbed = gsl_complex_rect(std::sqrt(1)+x, 0);
  

  // Create the Quantum Graph object
  QuantumGraphObject QGO;

  
  //Left Graph
  QGO.AddNode(right_circulator3); // 0
  QGO.AddNode(neumann4);    // 1
  QGO.AddNode(neumann4);    // 2
  QGO.AddNode(neumann4);    // 3
  QGO.AddNode(neumann3);    // 4
  
  //Right Graph
  QGO.AddNode(left_circulator3); // 5
  QGO.AddNode(neumann4);    // 6
  QGO.AddNode(neumann4);    // 7
  QGO.AddNode(neumann4);    // 8
  QGO.AddNode(neumann3);    // 9
  
  //Connection , no phase shift
  QGO.AddNode(left_circulator3); // 10
  QGO.AddNode(left_circulator3); // 11
  QGO.AddNode(neumann_reflector); // 12
  QGO.AddNode(neumann_reflector); // 13
  
  //Connection, pi phase shift
  QGO.AddNode(left_circulator3); // 15
  QGO.AddNode(left_circulator3); // 14
  QGO.AddNode(dirichlet_reflector); // 16
  QGO.AddNode(dirichlet_reflector); // 17
  
  
  // The Left and Right Graphs, have equal length bonds
  QGO.Connect(0, 1, length1);
  QGO.Connect(5, 6, length1_perturbed);
  
  QGO.Connect(0, 2, length2);
  QGO.Connect(5, 7, length2);
  
  QGO.Connect(0, 3, length3);
  QGO.Connect(5, 8, length3);
  
  QGO.Connect(1, 2, length4);
  QGO.Connect(6, 7, length4);
  
  QGO.Connect(2, 3, length5);
  QGO.Connect(7, 8, length5);
  
  QGO.Connect(1, 4, length6);
  QGO.Connect(6, 9, length6);
  
  QGO.Connect(2, 4, length7);
  QGO.Connect(7, 9, length7);
  
  QGO.Connect(3, 4, length8);
  QGO.Connect(8, 9, length8);
  
  //The connection between the graphs and the circulators
  QGO.Connect(1, 10, length9);
  QGO.Connect(8, 11, length9);
  
  QGO.Connect(3, 15, length9);
  QGO.Connect(6, 14, length9);
  
  //The short wire between the circulators.
  QGO.Connect(10, 11, length11);
  QGO.Connect(14, 15, length11);
  
  //The top circulators and the neumann reflectors
  QGO.Connect(10, 12, length10);
  QGO.Connect(11, 13, length10);
  
  //The bottom circulators and the dirichlet reflectors
  QGO.Connect(14, 16, length10);
  QGO.Connect(15, 17, length10);
  
  

  QGO.UpdateQuantumGraph();
  
  
  QGO.Normalize();


  return QGO;
 
}








