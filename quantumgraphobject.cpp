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
  nodes_[a]->AttachToBond(undirected_bonds_.back());
  nodes_[b]->AttachToBond(undirected_bonds_.back());
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
int QuantumGraphObject::GetIndexForBond(QuantumGraphBond* QGBP)
{
  // Take every undirected bond in our list.
  for (unsigned int ub=0; ub<undirected_bonds_.size(); ub++)
  {
    // First check to see if the forward bond has the pointer we're 
    // looking for.
    if ((undirected_bonds_[ub]->forward_bond()) == QGBP)
    {
      return 2 * ub;
    } // Then check the backward bond.
    else if ((undirected_bonds_[ub]->backward_bond()) == QGBP)
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
int QuantumGraphObject::GetIndexForNode(QuantumGraphNode* QGNP)
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
    QuantumGraphBond* forward_bond
        = undirected_bonds_[b]->forward_bond();

    // Our new Length vector contains the length information stored in
    // the directed bond we're referring to.
    gsl_vector_complex_set(newLVector, i, forward_bond->complex_length());


    // Our directed bond can tell us which node it starts at.
    QuantumGraphNode* FStartNode = forward_bond->start_node();

    // Our starting node can tell us which directed bonds arrive there.
    std::vector<QuantumGraphBond*> FIn 
        = FStartNode->GetIncomingBonds();

    // And we need to know what row in the node scattering matrix
    // corresponds to the outgoing bond.
    l = FStartNode->GetIndexForBond(forward_bond);

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
       int j = GetIndexForBond(FIn[m]);

       // Now actually build our graph scattering matrix's row i.
       gsl_complex Sij 
           = FStartNode->node_scattering_matrix_element(l, m);
       gsl_matrix_complex_set(newSMatrix, i, j, Sij);
    }

    // Now do exactly the same for the Backward Bond
    i = (2*b)+1; // Only difference: the backward bond is the next row
    QuantumGraphBond* backward_bond 
        = undirected_bonds_[b]->backward_bond();
    gsl_vector_complex_set(newLVector,i, 
                           backward_bond->complex_length());
    QuantumGraphNode* BStartNode = backward_bond->start_node();
    std::vector<QuantumGraphBond*> BIn 
        = BStartNode->GetIncomingBonds();
    l = BStartNode->GetIndexForBond(backward_bond);
    for (unsigned int m=0; m<BIn.size(); m++)
    {
       int j = GetIndexForBond(BIn[m]);
       gsl_complex Sij 
           = BStartNode->node_scattering_matrix_element(l, m);
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

// Does no checking for good input.
// Basically just copies whatever scattering matrix you gave it. It's up
// to the user to make sure this is set in a way that makes sense to the
// physical problem.
QuantumGraphNode::QuantumGraphNode(
    gsl_matrix_complex* node_scattering_matrix
)
{
  valence_ = node_scattering_matrix->size1;
  node_scattering_matrix_ 
      = gsl_matrix_complex_calloc(node_scattering_matrix->size1,
                                  node_scattering_matrix->size2);
  gsl_matrix_complex_memcpy(node_scattering_matrix_,
                            node_scattering_matrix);
}

// The default constructor means that there's nothing to describe and
// no connections to worry about
QuantumGraphNode::QuantumGraphNode()
{
  valence_ = 0;
  node_scattering_matrix_ = gsl_matrix_complex_calloc(0, 0);
}

// Have to delete the gsl_matrix_complex memory we allocated in the
// constructor.
QuantumGraphNode::~QuantumGraphNode()
{
  gsl_matrix_complex_free(node_scattering_matrix_);
}

// Just keep track of the pointers to the bonds the node's attached to 
// since the nodes are themselves children of the QuantumGraphObject.
void QuantumGraphNode::AttachToBond(QuantumGraphUndirectedBond* QGUB)
{
  attached_undirected_bonds_.push_back(QGUB);
}

// We don't keep track of the bonds that are attached to a node, only
// the undirected bonds. So we have to ask each undirected bond which
// of its bonds stores a pointer to this node as it's end_node_ member.
std::vector<QuantumGraphBond*> QuantumGraphNode::GetIncomingBonds()
{
  std::vector<QuantumGraphBond*> IncomingBonds;

  // Each undirected bond is attached to the node, so either the forward
  // or the backward bond ends at this node
  for (unsigned int m=0; m<attached_undirected_bonds_.size(); m++)
  {
    // If a bond ends at this node, it stores a pointer to this node.
    // by checking the end_node_ pointer for both the forward...
    if (attached_undirected_bonds_[m]->forward_bond()->end_node()
            == this)
    {
      IncomingBonds.push_back(attached_undirected_bonds_[m]
                                  ->forward_bond());
    } // ...and backward...
    else if (attached_undirected_bonds_[m]->backward_bond()->end_node()
                 == this)
    {
      IncomingBonds.push_back(attached_undirected_bonds_[m]
                                  ->backward_bond());
    }
    // ...bonds, we can make a vector of pointers to the incoming bonds
    // and return it to functions that need that information to build
    // the scattering matrix.
  }

  return IncomingBonds;
}


// Pretty straight forward. Iterate over all the attached undirected 
// bonds and see if the bond is in one of them.
int QuantumGraphNode::GetIndexForBond(QuantumGraphBond* QGBP)
{
  for (unsigned int m=0; m<attached_undirected_bonds_.size(); m++)
  {
    if (attached_undirected_bonds_[m]->HasDirectedBond(QGBP))
    {
      // Since the index of the bond in the scattering matrix is just
      // the position in the vector, return that position
      return m;
    }
  }
  
  // Failure returns something we can't use.
  return -1;
}


// Accessor Method... sort of... since I don't want to pass a pointer to
// data, which would remove encapsulation.
gsl_complex QuantumGraphNode::node_scattering_matrix_element(int l, int m)
{
  return gsl_matrix_complex_get(node_scattering_matrix_, l, m);
}


// Accessor Method
std::vector<QuantumGraphUndirectedBond*>
QuantumGraphNode::attached_undirected_bonds()
{
  return attached_undirected_bonds_;
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphBond Class
///
///

// Straightforward. Remember that we only store pointers that are 
// children of the QuantumGraphObject.
QuantumGraphBond::QuantumGraphBond(QuantumGraphNode* start_node, 
                                   QuantumGraphNode* end_node, 
                                   gsl_complex complex_length)
{
  start_node_ = start_node;
  end_node_ = end_node;
  complex_length_ = complex_length;
}

// Accessor Method
QuantumGraphNode* QuantumGraphBond::start_node()
{
  return start_node_;
}

// Accessor Method
QuantumGraphNode* QuantumGraphBond::end_node()
{
  return end_node_;
}

// Accessor Method
gsl_complex QuantumGraphBond::complex_length()
{
  return complex_length_;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphUndirectedBond Class
///
/// 


// Note that the undirected bond is the parent of two bonds. Other
// classes have to call this object's methods to get information about
// which bonds are stored here.
QuantumGraphUndirectedBond::QuantumGraphUndirectedBond(
    QuantumGraphNode* node_one,
    QuantumGraphNode* node_two, 
    gsl_complex complex_length
) :forward_bond_(node_one, node_two, complex_length), 
   backward_bond_(node_two, node_one, complex_length)
{
}


// 
bool QuantumGraphUndirectedBond::HasDirectedBond(QuantumGraphBond* QGB)
{
  if ((QGB == &forward_bond_) or (QGB == &backward_bond_))
  {
    return true;
  }
  else
  {
    return false;
  }
}


// Accessor method.
QuantumGraphBond* QuantumGraphUndirectedBond::forward_bond()
{
  return &forward_bond_;
}


// Accessor method.
QuantumGraphBond* QuantumGraphUndirectedBond::backward_bond()
{
  return &backward_bond_;
}

