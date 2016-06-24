

#include "../quantumgraphobject.h"

#include <iostream>

#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


int main()
{

  double beta = 0.3;
  std::clog << "Set Lengths" << std::endl;
  
  double Ltot = (1 + std::sqrt(2) + std::sqrt(3) 
                  + std::sqrt(5) + std::sqrt(6))/M_PI;
  
  gsl_complex length1 = gsl_complex_rect(std::sqrt(1)/Ltot, 0);
  gsl_complex length1_pos = gsl_complex_mul_real(length1, 1 + beta);
  gsl_complex length1_neg = gsl_complex_mul_real(length1, 1 - beta);
  
  gsl_complex length2 = gsl_complex_rect(std::sqrt(2)/Ltot, 0);
  gsl_complex length2_pos = gsl_complex_mul_real(length2, 1 + beta);
  gsl_complex length2_neg = gsl_complex_mul_real(length2, 1 - beta);
  
  gsl_complex length3 = gsl_complex_rect(std::sqrt(3)/Ltot, 0);
  gsl_complex length3_pos = gsl_complex_mul_real(length3, 1 + beta);
  gsl_complex length3_neg = gsl_complex_mul_real(length3, 1 - beta);
  
  gsl_complex length4 = gsl_complex_rect(std::sqrt(5)/Ltot, 0);
  gsl_complex length4_pos = gsl_complex_mul_real(length4, 1 + beta);
  gsl_complex length4_neg = gsl_complex_mul_real(length4, 1 - beta);
  
  gsl_complex length5 = gsl_complex_rect(std::sqrt(6)/Ltot, 0);
  gsl_complex length5_pos = gsl_complex_mul_real(length5, 1 + beta);
  gsl_complex length5_neg = gsl_complex_mul_real(length5, 1 - beta);


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

  /// Define 1-bond Neumann node (open reflector)
  gsl_matrix_complex* neumann_reflector = gsl_matrix_complex_calloc(1,1);
  gsl_matrix_complex_set(neumann_reflector, 0, 0, GSL_COMPLEX_ONE);

  /// Define 1-bond Dirichlet node (closed reflector)
  gsl_matrix_complex* dirichlet_reflector = gsl_matrix_complex_calloc(1,1);
  gsl_matrix_complex_set(dirichlet_reflector, 0, 0, GSL_COMPLEX_NEGONE);



  std::clog << "Create Quantum Graph Object" << std::endl;
  QuantumGraphObject QGO;

  std::clog << "Add Nodes to Object" << std::endl;
  QGO.AddNode(neumann_reflector);
  QGO.AddNode(neumann4);
  QGO.AddNode(neumann3);
  QGO.AddNode(neumann4);
  QGO.AddNode(neumann3);
  QGO.AddNode(neumann4);
  QGO.AddNode(dirichlet_reflector);

  std::clog << "Connect Nodes with Bonds" << std::endl;
  QGO.Connect(1, 0, length1_pos, length1_neg);
  QGO.Connect(2, 1, length2_pos, length2_neg);
  QGO.Connect(3, 1, length3_pos, length3_neg);
  QGO.Connect(4, 1, length4_pos, length4_neg);
  QGO.Connect(3, 2, length5_pos, length5_neg);
  QGO.Connect(4, 3, length5_pos, length5_neg);
  QGO.Connect(5, 2, length4_pos, length4_neg);
  QGO.Connect(5, 3, length3_pos, length3_neg);
  QGO.Connect(5, 4, length2_pos, length2_neg);
  QGO.Connect(6, 5, length1_pos, length1_neg);
  


  std::clog << "Updating Quantum Graph Matrices" << std::endl;
  QGO.UpdateQuantumGraph();

  std::clog << "Printing the Matrix Information" << std::endl;
  std::cout << QGO << std::endl;


  return 0;
}
