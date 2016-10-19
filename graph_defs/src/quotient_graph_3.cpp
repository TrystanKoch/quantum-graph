

#include "quantumgraphobject.h"

#include <iostream>

#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


int main()
{

  double Ltot = 2*(1 + std::sqrt(2) + std::sqrt(3) 
                     + std::sqrt(5) + std::sqrt(6)
                     + 2*std::sqrt(7) + 2*std::sqrt(10)
                  )/M_PI;
  
  gsl_complex length1 = gsl_complex_rect(std::sqrt(1)/Ltot, 0);
  
  gsl_complex length2 = gsl_complex_rect(std::sqrt(2)/Ltot, 0);
  
  gsl_complex length3 = gsl_complex_rect(std::sqrt(3)/Ltot, 0);
  
  gsl_complex length4 = gsl_complex_rect(std::sqrt(5)/Ltot, 0);
  
  gsl_complex length5 = gsl_complex_rect(std::sqrt(6)/Ltot, 0);
  
  gsl_complex length6 = gsl_complex_rect(std::sqrt(7)/Ltot, 0);
  
  gsl_complex length7 = gsl_complex_rect(std::sqrt(10)/Ltot, 0);


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
  
  /// Define 4-bond Circulator node
  // Looks weird, because rows 1 and 2 are across from each other.
  // Scattering is from 1->3,3->2,2->4,4->1
  std::vector<std::vector<double>> C4 = 
  {
    { 0,  0,  0,  1},
    { 0,  0,  1,  0},
    { 1,  0,  0,  0},
    { 0,  1,  0,  0}
  };

  gsl_matrix_complex* circulator4 = gsl_matrix_complex_calloc(4, 4);
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
    { 
      gsl_matrix_complex_set(circulator4, i, j, gsl_complex_rect(C4[i][j], 0));
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
  QGO.AddNode(neumann3);    // 0
  QGO.AddNode(neumann3);    // 1
  QGO.AddNode(left_circulator3); // 2
  QGO.AddNode(neumann3);    // 3
  
  //Right Graph
  QGO.AddNode(neumann3);    // 4
  QGO.AddNode(neumann3);    // 5
  QGO.AddNode(right_circulator3); // 6
  QGO.AddNode(neumann3);    // 7
  
  //Connection , no phase shift
  QGO.AddNode(neumann_reflector); // 8
  QGO.AddNode(circulator4);       // 9
  QGO.AddNode(neumann_reflector); // 10
  
  //Connection, pi phase shift
  QGO.AddNode(dirichlet_reflector); // 11
  QGO.AddNode(circulator4);         // 12
  QGO.AddNode(dirichlet_reflector); // 13
  

  std::clog << "Connect Nodes with Bonds" << std::endl;
  
  // The Left and Right Graphs, have equal length bonds
  QGO.Connect(0, 1, length1);
  QGO.Connect(4, 5, length1);
  
  QGO.Connect(0, 2, length2);
  QGO.Connect(4, 6, length2);
  
  QGO.Connect(0, 3, length3);
  QGO.Connect(4, 7, length3);
  
  QGO.Connect(1, 2, length4);
  QGO.Connect(5, 6, length4);
  
  QGO.Connect(2, 3, length5);
  QGO.Connect(6, 7, length5);
  
  //The connection between the graphs and the circulators
  QGO.Connect(1, 9, length6);
  QGO.Connect(7, 9, length6);
  
  QGO.Connect(3, 12, length6);
  QGO.Connect(5, 12, length6);
  
  //The top circulator and the neumann reflectors
  QGO.Connect(9,  8, length7);
  QGO.Connect(9, 10, length7);
  
  //The bottom circulator and the dirichlet reflectors
  QGO.Connect(12, 11, length7);
  QGO.Connect(12, 13, length7);
  


  std::clog << "Updating Quantum Graph Matrices" << std::endl;
  QGO.UpdateQuantumGraph();

  std::clog << "Printing the Matrix Information" << std::endl;
  std::cout << QGO << std::endl;

  std::clog << Ltot;

  return 0;
}
