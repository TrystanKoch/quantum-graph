/*----------------------------------------------------------------------
   quantumgraphobject.cpp 
     - Function declarations for a class representing a quantum graph 
       whose structure resembles the quantum graph
     - 

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
//    QuantumGraphObject
//
//  QuantumGraphObject
//    represents a Quantum Graph's Structure

#ifndef QUANTUMGRAPHOBJECT_H_
#define QUANTUMGRAPHOBJECT_H_

#include "quantumgraph.h"

#include <ostream>
#include <vector>
#include <string>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class QuantumGraphObject;
class QuantumGraphNode;
class QuantumGraphUndirectedBond;
class QuantumGraphBond;

// 
//
//
class QuantumGraphObject: public QuantumGraph
{
  private:
    std::vector<QuantumGraphNode*> Nodes;
    std::vector<QuantumGraphUndirectedBond*> UndirectedBonds;
  public:
    QuantumGraphObject();
    QuantumGraphObject(const QuantumGraphObject&);
    ~QuantumGraphObject();
    QuantumGraphObject(std::vector<QuantumGraphUndirectedBond*>, 
                       std::vector<QuantumGraphNode*>);
    void AddNode();
    void AddNode(gsl_matrix_complex*);
    void Connect(QuantumGraphNode&, QuantumGraphNode&, gsl_complex);
    void Connect(unsigned int, unsigned int, gsl_complex);
    void UpdateQuantumGraph();
    int GetBondIndexFromPointer(QuantumGraphBond*);
    int GetNodeIndexFromPointer(QuantumGraphNode*);
};



class QuantumGraphNode
{
  private:
    std::vector<QuantumGraphUndirectedBond*> ConnectedBonds;
    gsl_matrix_complex* SMatrix;
    unsigned int Valence;
  public:
    QuantumGraphNode();
    QuantumGraphNode(gsl_matrix_complex*);
    ~QuantumGraphNode();
    //QuantumGraphNode(const QuantumGraphNode&);
    //std::vector<QuantumGraphNode&> GetConnectedNodes() const;
    //std::vector<QuantumGraphBond&> GetOutgoingBonds() const;
    std::vector<QuantumGraphBond*> GetIncomingBonds();
    std::vector<QuantumGraphUndirectedBond*> GetConnectedUBonds();
    int GetBondIndexAtNodeFromPointer(QuantumGraphBond*);
    void ConnectToBond(QuantumGraphUndirectedBond*);
    gsl_complex GetMatrixElement(int, int);
};



class QuantumGraphBond
{
  private:
    gsl_complex ComplexLength;
    QuantumGraphNode* StartNode;
    QuantumGraphNode* EndNode;
  public:
    //QuantumGraphBond();
    //QuantumGraphBond(const QuantumGraphBond&);
    QuantumGraphBond(QuantumGraphNode*, QuantumGraphNode*, gsl_complex);
    QuantumGraphNode* GetStartNode();
    QuantumGraphNode* GetEndNode();
    gsl_complex GetBondLength();
};


class QuantumGraphUndirectedBond
{
  private:
    QuantumGraphBond ForwardBond;
    QuantumGraphBond BackwardBond;
    void CheckBondLengths();
  public:
    //QuantumGraphUndirectedBond();
    //QuantumGraphUndirectedBond(const QuantumGraphUndirectedBond&);
    QuantumGraphUndirectedBond(QuantumGraphNode*, QuantumGraphNode*, gsl_complex);
    QuantumGraphBond* GetForwardBond();
    QuantumGraphBond* GetBackwardBond();
    bool HasDirectedBond(QuantumGraphBond*);
    //std::vector<QuantumGraphNode&> GetAttachedNodes const;
    
};







#endif
