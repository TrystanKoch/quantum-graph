

#include "quantumgraphobject.h"

#include <iostream>

#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


int main()
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

  gsl_complex length10_per = gsl_complex_rect(std::sqrt(14)+10, 0);


  std::clog << "Set Node Scattering Matrices" << std::endl;
  
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



  std::clog << "Create Quantum Graph Object" << std::endl;
  QuantumGraphObject QGO;

  std::clog << "Add Nodes to Object" << std::endl;
  
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
  

  std::clog << "Connect Nodes with Bonds" << std::endl;
  
  // The Left and Right Graphs, have equal length bonds
  QGO.Connect(0, 1, length1);
  QGO.Connect(5, 6, length1);
  
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
  QGO.Connect(11, 13, length10_per);
  
  //The bottom circulators and the dirichlet reflectors
  QGO.Connect(14, 16, length10);
  QGO.Connect(15, 17, length10);
  


  std::clog << "Updating Quantum Graph Matrices" << std::endl;
  QGO.UpdateQuantumGraph();
  
  std::clog << "Normalizing bond lengths" <<std::endl;
  QGO.Normalize();

  std::clog << "Printing the Matrix Information" << std::endl;
  std::cout << QGO << std::endl;
  
  return 0;
}
