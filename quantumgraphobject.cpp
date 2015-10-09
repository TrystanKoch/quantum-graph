/*----------------------------------------------------------------------
   quantumgraphobject.cpp 
     - Function definitions for a class representing the structure of a
       quantum graph.
     - Update function for the quantum graph class. A protected method
       of the QuantumGraph class that creates a useable mathematical
       representation of the graph we are modeling.
     - Function definitions for classes modeling nodes_, bonds, and 
       undirected bonds in the quantum graph's structure

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

//  Provides source code for the classes
//    QuantumGraphObject
//    QuantumGraphNode
//    QuantumGraphBond
//    QuantumGraphUndirectedBond
//
//  You can find more information in the header file. The comments here
//  regard the implementation itself, and how the methods work.

#include "quantumgraph.h"
#include "quantumgraphobject.h"

#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphObject Class
///
/// 


// Without input, set both lists to empty.
QuantumGraphObject::QuantumGraphObject()
{
  nodes_ = {};
  undirected_bonds_ = {};
}


// Dynamically add a new node. Not necessarily connected to anything.
// The scattering matrix supplied could be anything, as there is no
// check on what gets input.
void QuantumGraphObject::AddNode(gsl_matrix_complex* S)
{
  nodes_.push_back(new QuantumGraphNode(S));
}


// Connect node a with node b, using an undirected bond of complex 
// length L. Can only add bonds this way, so bonds are always connected.
void QuantumGraphObject::Connect(unsigned int a, unsigned int b, 
                                 gsl_complex L)
{
  undirected_bonds_.push_back(
               new QuantumGraphUndirectedBond(nodes_[a], nodes_[b], L));
  nodes_[a]->ConnectToBond(undirected_bonds_.back());
  nodes_[b]->ConnectToBond(undirected_bonds_.back());
}



// Loops through both member vectors and ensures that the memory 
// allocated to the stored pointer is deleted. Important because we're
// dynamically adding nodes and bonds.
QuantumGraphObject::~QuantumGraphObject()
{
  // Free nodes from memory
  for (unsigned int i=0; i<nodes_.size(); i++)
  {
    delete nodes_[i];
  }

  // Free bonds from memory. QuantumGraphUndirectedBond class should
  // handle the directed bonds itself.
  for (unsigned int i=0; i<undirected_bonds_.size(); i++)
  {
    delete undirected_bonds_[i];
  }
}


// We only store the bonds as a list, so if we find a pointer to one
// somewhere else -- say, a node asking for all the incoming bonds --
// we need to know what the index of that bond is.
//
// Note that in this convention, the index "i" of the forward bond in
// the undirected bond "b" is i=2*b, while the backward bond is given
// index i=2*b + 1. 
int QuantumGraphObject::GetBondIndexFromPointer(QuantumGraphBond* QGBP)
{
  // Take every undirected bond in our list.
  for (unsigned int ub=0; ub<undirected_bonds_.size(); ub++)
  {
    // First check to see if the forward bond has the pointer we're 
    // looking for.
    if ((undirected_bonds_[ub]->GetForwardBond()) == QGBP)
    {
      return 2 * ub;
    } // Then check the backward bond.
    else if ((undirected_bonds_[ub]->GetBackwardBond()) == QGBP)
    {
      return 2 * ub + 1;
    }
  }

  // Failure gives us a negative index.
  return -1;
}


// We only store the nodes as a list, so if we find a pointer to one
// somewhere else we need to know what the index of that bond is.
//
// less important than the bond version, but still useful.
int QuantumGraphObject::GetNodeIndexFromPointer(QuantumGraphNode* QGNP)
{
  for (unsigned int n=0; n<nodes_.size(); n++)
  {
    if (nodes_[n] == QGNP)
    {
      return n;
    }
  }

  // Failure gives us a negative index.
  return -1;
}



// Once we're done making the structure, we have to pass it back to the
// base class data members before we can use it in modeling.
void QuantumGraphObject::UpdateQuantumGraph()
{
  // This matrix will hold the scattering matrix as we are building it.
  // At the end of the file it will be the scattering matrix containing
  // bond-to-bond scattering for the whole graph. It has to have doubled
  // dimensions because it is a graph for directed bonds and we have the
  // number of undirected bonds.
  gsl_matrix_complex* newSMatrix 
      = gsl_matrix_complex_calloc(2 * undirected_bonds_.size(),
                                  2 * undirected_bonds_.size());

  // Similarly, this vector will be our new complex length vector. But
  // Before we pass it to the protected setter function, we have to
  // build it from the information we have about the structure of the
  // quantum graph.
  gsl_vector_complex* newLVector
      = gsl_vector_complex_calloc(2 * undirected_bonds_.size());


  // For every undirected bond in the graph, determine which node each
  // directed bond starts at. Find which directed bonds end at that
  // node. Then determine which element in the node's scattering matrix
  // corresponds to the scattering from the latter to the former. 
  // Finally, store this in the graph's scattering matrix at the correct
  // position.
  // 
  // Harder than it seems: only the QuantumGraphObject knows which
  // matrix columns/rows correspond to which bond, and only through its
  // position in the vector. Similarly, only the bonds know which nodes
  // they are connected to, and only the node knows which other bonds
  // arrive at the node and which column/row of the node's scattering
  // matrix that the bond corresponds to. Since the node scattering
  // matrix is a different size than the graph's scattering matrix, we
  // have to sort out the right indices.
  //
  // I use the following:
  //  b - The undirected bond's position in the list
  //  i - The row of the graph scatttering matrix corresponding to the
  //      appropriate directed bond of undirected bond "b".
  //  j - The column of the graph scattering matrix corresponding to the
  //      incoming directed bond, scattering into the appropriate 
  //      directed bond of undirected bond "b"
  //  l - Row of the Node scattering matrix corresponding to "i"
  //  m - Column of the Node scattering matrix corresponding to "j"
  //
  for (unsigned int b=0; b<undirected_bonds_.size(); b++)
  {
    // These will correspond to the directed bond we're scattering to.
    // We build the matrix row by row.
    int i, l;

    // First consider the forward bond of b.
    // It will have index 2*b in the graph scattering matrix.
    i = 2*b; 

    // We need to know what directed bond row i refers to. We get the
    // pointer that refers to the object.
    QuantumGraphBond* ForwardBond 
        = undirected_bonds_[b]->GetForwardBond();

    // Our new Length vector contains the length information stored in
    // the directed bond we're referring to.
    gsl_vector_complex_set(newLVector, i, ForwardBond->GetBondLength());


    // Our directed bond can tell us which node it starts at.
    QuantumGraphNode* FStartNodeP = ForwardBond->GetStartNode();

    // Our starting node can tell us which directed bonds arrive there.
    std::vector<QuantumGraphBond*> FIn 
        = FStartNodeP->GetIncomingBonds();

    // And we need to know what row in the node scattering matrix
    // corresponds to the outgoing bond.
    l = FStartNodeP->GetBondIndexAtNodeFromPointer(ForwardBond);

    // For each bond incoming to the start node, we find its 
    // corresponding graph scattering matrix index j. Luckily, the order
    // of the incoming bonds stored in the node is already the index
    // of the incoming bond, so we can loop over the column index m.
    for (unsigned int m=0; m<FIn.size(); m++)
    {
       // Figure out where the bond is in the QuantumGraphObject list.
       // Whether the "forward" or "backward" vector is considered
       // before or after the other in this list does not matter as long
       // as we're consistent. For now, I assume the Forward bond is the
       // first one.
       int j = GetBondIndexFromPointer(FIn[m]);

       // Now actually build our graph scattering matrix's row i.
       gsl_complex Sij = FStartNodeP->GetMatrixElement(l, m);
       gsl_matrix_complex_set(newSMatrix, i, j, Sij);
    }

    // Now Do exactly the same for the Backward Bond
    i = (2*b)+1; // Only difference: the backward bond is the next row
    QuantumGraphBond* BackwardBond 
        = undirected_bonds_[b]->GetBackwardBond();
    gsl_vector_complex_set(newLVector,i,BackwardBond->GetBondLength());
    QuantumGraphNode* BStartNodeP = BackwardBond->GetStartNode();
    std::vector<QuantumGraphBond*> BIn 
        = BStartNodeP->GetIncomingBonds();
    l = BStartNodeP->GetBondIndexAtNodeFromPointer(BackwardBond);
    for (unsigned int m=0; m<BIn.size(); m++)
    {
       int j = GetBondIndexFromPointer(BIn[m]);
       gsl_complex Sij = BStartNodeP->GetMatrixElement(l, m);
       gsl_matrix_complex_set(newSMatrix, i, j, Sij);
    }

    // Probably possible to have a helper function to clean this up,
    // especially since much of the same work is done for both forward
    // and backward bond. 

  } 
  // Got through all the undirected bonds in the list.

  // Protected function of the QuantumGraph class. But we're in a 
  // derived class so we can use it.
  setGraph(newLVector, newSMatrix);

  // Clean the workspace of our temporary matrix and vector.
  gsl_matrix_complex_free(newSMatrix);
  gsl_vector_complex_free(newLVector);
}







////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphNode Class
///
///   
QuantumGraphNode::QuantumGraphNode(gsl_matrix_complex* SS)
{
  Valence = SS->size1;
  SMatrix = gsl_matrix_complex_calloc(SS->size1, SS->size2);
  gsl_matrix_complex_memcpy(SMatrix, SS);
}

QuantumGraphNode::QuantumGraphNode()
{
  Valence = 0;
  SMatrix = gsl_matrix_complex_calloc(0, 0);
}

QuantumGraphNode::~QuantumGraphNode()
{
  gsl_matrix_complex_free(SMatrix);
}

void QuantumGraphNode::ConnectToBond(QuantumGraphUndirectedBond* QGUB)
{
  ConnectedBonds.push_back(QGUB);
}


std::vector<QuantumGraphBond*> QuantumGraphNode::GetIncomingBonds()
{
  std::vector<QuantumGraphBond*> InBondsPVec;
  for (unsigned int m=0; m<ConnectedBonds.size(); m++)
  {
    if (ConnectedBonds[m]->GetForwardBond()->GetEndNode() == this)
    {
      InBondsPVec.push_back(ConnectedBonds[m]->GetForwardBond());
    }
    else if (ConnectedBonds[m]->GetBackwardBond()->GetEndNode() == this)
    {
      InBondsPVec.push_back(ConnectedBonds[m]->GetBackwardBond());
    }
  }

  return InBondsPVec;
}


int 
QuantumGraphNode::GetBondIndexAtNodeFromPointer(QuantumGraphBond* QGBP)
{
  for (unsigned int m=0; m<ConnectedBonds.size(); m++)
  {
    if (ConnectedBonds[m]->HasDirectedBond(QGBP))
    {
      return m;
    }
  }
  
  return -10;
}



gsl_complex QuantumGraphNode::GetMatrixElement(int l, int m)
{
  return gsl_matrix_complex_get(SMatrix, l, m);
}

std::vector<QuantumGraphUndirectedBond*>
QuantumGraphNode::GetConnectedUBonds()
{
  return ConnectedBonds;
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphBond Class
///
/// 
QuantumGraphBond::QuantumGraphBond(QuantumGraphNode* Start, 
                                   QuantumGraphNode* End, 
                                   gsl_complex LL)
{
  StartNode = Start;
  EndNode = End;
  ComplexLength = LL;
}

QuantumGraphNode* QuantumGraphBond::GetStartNode()
{
  return StartNode;
}

QuantumGraphNode* QuantumGraphBond::GetEndNode()
{
  return EndNode;
}

gsl_complex QuantumGraphBond::GetBondLength()
{
  return ComplexLength;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphUndirectedBond Class
///
///  
QuantumGraphUndirectedBond::QuantumGraphUndirectedBond(
  QuantumGraphNode* ForwardStart,
  QuantumGraphNode* ForwardEnd, 
  gsl_complex LL
) :ForwardBond(ForwardStart, ForwardEnd, LL), 
   BackwardBond(ForwardEnd, ForwardStart, LL)
{
}



bool QuantumGraphUndirectedBond::HasDirectedBond(QuantumGraphBond* QGB)
{
  if ((QGB == &ForwardBond) or (QGB == &BackwardBond))
  {
    return true;
  }
  else
  {
    return false;
  }
}

QuantumGraphBond* QuantumGraphUndirectedBond::GetForwardBond()
{
  return &ForwardBond;
}

QuantumGraphBond* QuantumGraphUndirectedBond::GetBackwardBond()
{
  return &BackwardBond;
}

