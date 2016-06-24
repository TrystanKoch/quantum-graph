

#include "../quantumgraphobject.h"

#include <iostream>

#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

const unsigned int B = 12;

const std::vector<std::vector<double>> SS = {
  {0, 0, 0, 0, 0, 0, -1./3., 2./3., 2./3., 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 2./3., -1./3., 2./3., 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 2./3., 2./3., -1./3., 0, 0, 0},
  {0, 0, 2./3., 0, 0, 2./3., 0, 0, 0, -1./3., 0, 0},
  {2./3., 0, 0, 2./3., 0, 0, 0, 0, 0, 0, -1./3., 0},
  {0, 2./3., 0, 0, 2./3., 0, 0, 0, 0, 0, 0, -1./3.},
  {-1./3., 0, 0, 2./3., 0, 0, 0, 0, 0, 0, 2./3., 0},
  {0, -1./3., 0, 0, 2./3., 0, 0, 0, 0, 0, 0, 2./3.},
  {0, 0, -1./3., 0, 0, 2./3., 0, 0, 0, 2./3., 0, 0},
  {2./3., 0, 0, -1./3., 0, 0, 0, 0, 0, 0, 2./3., 0},
  {0, 2./3., 0, 0, -1./3., 0, 0, 0, 0, 0, 0, 2./3.},
  {0, 0, 2./3., 0, 0, -1./3., 0, 0, 0, 2./3., 0, 0}
};


//Side lengths of Tetrahedron with positive volume
const std::vector<double> LL = {
  1. + M_1_PI,
  sqrt(5) + M_1_PI,
  sqrt(3) + M_1_PI,
  sqrt(2) + M_1_PI,
  sqrt(6) + M_1_PI,
  sqrt(7) + M_1_PI,
  1. - M_1_PI,
  sqrt(5) - M_1_PI,
  sqrt(3) - M_1_PI,
  sqrt(2) - M_1_PI,
  sqrt(6) - M_1_PI,
  sqrt(7) - M_1_PI
};


int main()
{


  gsl_matrix_complex* S = gsl_matrix_complex_calloc(12, 12);
  for (int i=0; i<12; i++)
  {
    for (int j=0; j<12; j++)
    { 
      gsl_matrix_complex_set(S, i, j, gsl_complex_rect(SS[i][j], 0));
    }
  }

  gsl_vector_complex* L = gsl_vector_complex_calloc(12);
  for (int i=0; i<12; i++)
  {
    gsl_vector_complex_set(L, i, gsl_complex_rect(LL[i], 0));
  }

  std::clog << "Create Quantum Graph Object" << std::endl;
  QuantumGraph QG(L,S);


  std::clog << "Printing the Matrix Information" << std::endl;
  std::cout << QG << std::endl;


  return 0;
}
