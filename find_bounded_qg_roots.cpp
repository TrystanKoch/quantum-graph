/*----------------------------------------------------------------------
   find_bounded_qg_roots.cpp 
     - A program to take boundaries that have been found to contain only
       one root and quickly find the root to the best possible accuracy
       using Newton's Method.

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

#include "quantumgraph.h"
#include "quantumgraphrootfinding.h"

#include <sstream>
#include <limits>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_qrng.h>


int main(int argc, char* argv[])
{

  std::vector<std::string> header;
  std::string line;
  while (std::getline(std::cin, line))
  {
    if (line[0]=='#')
    {
      header.push_back(line);
      if (!line.compare("# END"))
      {
        break;
      }
    }
  }

  QuantumGraph QG = QuantumGraph(header);

  std::stringstream ss;

  BoundingBox BB;
  gsl_complex root;
  gsl_qrng* qrng = gsl_qrng_alloc(gsl_qrng_sobol, 2);
  
  std::cout.precision(std::numeric_limits<double>::max_digits10);
  std::cout << std::scientific << std::showpos;

  unsigned long long int count = 0;
  while (std::getline(std::cin, line))
  {
    ss.clear();
    ss.str(line);
    ss >> BB.realMin;
    ss >> BB.realMax;
    ss >> BB.imagMin;
    ss >> BB.imagMax;
    
    root = newtonRootFinderQRNG(BB, qrng, QG);
    std::cout << GSL_REAL(root) << "  " << GSL_IMAG(root) << "    " << gsl_complex_abs(QG.characteristic(root))<<std::endl;
    count++;
  }
  std::clog << "Roots Polished: " << count << std::endl;
  return 0;

}
