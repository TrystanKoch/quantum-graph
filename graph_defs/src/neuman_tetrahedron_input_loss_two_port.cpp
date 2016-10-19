

#include "quantumgraphobject.h"

#include <iostream>

#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


int main()
{

  std::clog << "Set Lengths" << std::endl;
  gsl_complex length1 = gsl_complex_rect(std::sqrt(1), 0);
  gsl_complex length2 = gsl_complex_rect(std::sqrt(5), 0);
  gsl_complex length3 = gsl_complex_rect(std::sqrt(3), 0);
  gsl_complex length4 = gsl_complex_rect(std::sqrt(2), 0);
  gsl_complex length5 = gsl_complex_rect(std::sqrt(6), 0);
  gsl_complex length6 = gsl_complex_rect(std::sqrt(7), 0);

  std::clog << "Set Node Scattering Matrix" << std::endl;
  std::vector<std::vector<double>> N3 = 
  {
    {-1./3.,  2./3.,  2./3.},
    { 2./3., -1./3.,  2./3.},
    { 2./3.,  2./3., -1./3.}
  };

  std::vector<std::vector<double>> N3_loss = 
  {
    {-1./2.,  1./2.,  1./2.},
    { 1./2., -1./2.,  1./2.},
    { 1./2.,  1./2., -1./2.}
  };

  gsl_matrix_complex* Neumann3 = gsl_matrix_complex_calloc(3, 3);
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    { 
      gsl_matrix_complex_set(Neumann3, i, j, 
                             gsl_complex_rect(N3[i][j], 0));
    }
  }

  gsl_matrix_complex* Neumann3_loss = gsl_matrix_complex_calloc(3, 3);
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    { 
      gsl_matrix_complex_set(Neumann3_loss, i, j, 
                             gsl_complex_rect(N3_loss[i][j], 0));
    }
  }

  std::clog << "Create Quantum Graph Object" << std::endl;
  QuantumGraphObject QGO;

  std::clog << "Add Nodes to Object" << std::endl;
  QGO.AddNode(Neumann3_loss);
  QGO.AddNode(Neumann3_loss);
  QGO.AddNode(Neumann3);
  QGO.AddNode(Neumann3);

  std::clog << "Connect Nodes with Bonds" << std::endl;
  QGO.Connect(0, 1, length1);
  QGO.Connect(0, 2, length2);
  QGO.Connect(0, 3, length3);
  QGO.Connect(3, 1, length4);
  QGO.Connect(1, 2, length5);
  QGO.Connect(2, 3, length6);

  std::clog << "Updating Quantum Graph Matrices" << std::endl;
  QGO.UpdateQuantumGraph();

  std::clog << "Printing the Matrix Information" << std::endl;
  std::cout << QGO << std::endl;


  return 0;
}
