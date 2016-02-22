/*----------------------------------------------------------------------
   openquantumgraphobject.h
     - Class definition and Method prototypes for a class representing 
       the structure of a quantum graph.
     - Class Definitions and Function prototypes for classes modeling
       nodes, bonds, and undirected bonds in the quantum graph's 
       structure.

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

#ifndef OpenQuantumGraphObject_H_
#define OpenQuantumGraphObject_H_

class QuantumGraphLead
{
  private:
    
  public:
    
}

class QuantumGraphPort : public QuantumGraphNode
{
  private:
    QuantumGraphLead* attached_lead_;
    gsl_complex external_reflection_coefficient_;
    gsl_vector_complex* outgoing_coupling_coefficients_;
    gsl_vector_complex* incoming_coupling_coefficients_;
  public:
    QuantumGraphPort(const gsl_matrix_complex*);
    ConnectToLead(const unsigned int);
}
