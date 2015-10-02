
#ifndef TETRAHEDRON_GRAPH_H_
#define TETRAHEDRON_GRAPH_H_

#include <vector>

#include <cmath>

////////////////////////////////////////////////////////////////////////
/// Tetrahedron Definitions
//*

// Tetrahedron Scattering Matrix real part (no Loss)
const std::vector<std::vector<double>> SSreal = {
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

// Tetrahedron Scattering Matrix imaginary part (no Loss)
const std::vector<std::vector<double>> SSimag = {
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
};

// Side lengths of Tetrahedron with positive volume
const std::vector<double> LLreal = {
  1.,
  std::sqrt(5),
  std::sqrt(3),
  std::sqrt(2),
  std::sqrt(6),
  std::sqrt(7),
  1.,
  std::sqrt(5),
  std::sqrt(3),
  std::sqrt(2),
  std::sqrt(6),
  std::sqrt(7)
};

// Side lengths of Tetrahedron with positive volume (no Loss)
const std::vector<double> LLimag = {
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0
};

// Known root from previous program. Gives small vector norm.
//const gsl_complex ROOT1 = gsl_complex_rect(+8.85545367966650288e-01,0);
//const gsl_complex ROOT2 = gsl_complex_rect(+9.74748946419686457e-01,0);
//const gsl_complex ROOT3 = gsl_complex_rect(+1.09097781923122095e+00,0);
//const gsl_complex ROOT4 = gsl_complex_rect(+1.64642887492075896e+00,0);

#endif
