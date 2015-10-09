/*----------------------------------------------------------------------
   quantumgraphobject.h
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



//  Provides the classes
//    QuantumGraphObject
//    QuantumGraphNode
//    QuantumGraphBond
//    QuantumGraphUndirectedBond
//
//  QuantumGraphObject
//    represents a Quantum Graph's structure.
//    
//  QuantumGraphNode
//    
//    
//  QuantumGraphBond
//    
//    
//  QuantumGraphUndirectedBond
//    
//

//  For further documentation for the classes declared herein, there are
//  additional comments at the top of the class definitions.


//  Organizational Note:
//    I initially wished to seperate the header and the source code for
//    each of these four classes into eight seperate header/source files
//    but this cluttered the workspace a bit too much. Because these
//    classes are really only designed to work with each other, it makes
//    more sense to keep them in one source and one header file from
//    an organizational standpoint. This also simplifies compliation
//    commands drastically.
//
//    The downside to this is that this file and the source file are
//    both exceedingly long. I have attempted to make it clear where the
//    definitions and source code for one class ends and the next begins
//    but given that comments are comments, and these files contain many
//    comments, it might be hard to see these seperation points, which
//    I've marked with multiple full lines of forward slashes for the
//    best visibility. Hopefully, this does not excessively hinder your
//    exploration of the headers and source files.
//
//    To make things a bit neater, I have also tried to keep prototypes
//    and source code definitions for methods in the same order. Most 
//    documentation for use of these objects in programs is in this file
//    while specific comments on how the code operates can be found in
//    the source itself.


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// BEGINNING OF HEADER FILE

#ifndef QUANTUMGRAPHOBJECT_H_
#define QUANTUMGRAPHOBJECT_H_

#include "quantumgraph.h"

#include <ostream>
#include <vector>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// Forward Declarations
///
///   The three classes are intertwinned: we need forward declarations
///   so that they know when their vectors of pointers to the other 
///   classes are well defined.

class QuantumGraphObject;
class QuantumGraphNode;
class QuantumGraphUndirectedBond;
class QuantumGraphBond;







////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphObject Class
///
///
///  Is a QuantumGraph. 
///   The QuantumGraph class contains the mathematical representation 
///   and computational methods for our model. This QuantumGraphObject 
///   class serves primarily to create models in a more intuitive
///   manner, letting us specify nodes and connecting bonds instead of 
///   having to write out the full graph scattering matrix by hand for 
///   every new problem.
///
///
///  Conceptually models the structre of a Quantum Graph.
///   A graph is a collection of nodes, sometimes named vertices,
///   connected by bonds, sometimes named edges. A quantum graph is a
///   graph with a measure and a differential operator defined on the
///   bonds, and scattering conditions defined at the nodes.
///
///   The nodes contain the boundary conditions of the problem. Each
///   node is associated with a scattering matrix that relates a
///   function on a bond to the functions on bonds connected to it.
///   The Nodes are discussed further in the class QuantumGraphNode.
///
///   Each connection between nodes has a length. It may also have a 
///   constant loss or constant potential associated with it, and it
///   might even have a magnetic potential defined on it. These 
///   connections are termed "undirected bonds," and act much like we
///   expect a one-dimensional spatial domain might under physical 
///   interpretations of the differential operator. The Undirected Bonds
///   are discussed further in the class QuantumGraphUndirectedBond.
///
///   Because scattering problems are easiest to consider when we look 
///   at a directed graph, and any linear problem on an undirected graph
///   may be broken down into two directed bonds pointing in opposite
///   directions, we consider the basic bond structure of our Quantum
///   Graph Object not to be the undirected bond, but the directed bond,
///   which, for convenience, we simply term "bond". The Bonds are
///   discussed further in the class QuantumGraphBond.
///
///   By creating nodes with boundary conditions we can specify 
///   ourselves, based on the particular problem we wish to solve, then 
///   linking them together via undirected bonds we have the structure
///   of a quantum graph. The directed bonds are fundamental to the
///   model, but because they are always in opposing pairs, we consider
///   them as a sub unit of the undirected bonds when we declare a 
///   problem.
///
///   There are many possible types of bonds and nodes. It is expected 
///   that more specialized derived classes will simplify adding nodes
///   and bonds of verious types for various problems. These derived
///   classes can be directly used in the QuantumGraphObject class
///   because C++ does that intrinsically.
///
///
///  Examples of use:
///
///  -General-
///   Generally, this object is designed to specify a particular graph
///   to be modelled. Other programs are intended to use the output of
///   a program using this graph. Though it is quite possible to both
///   create the graph's scattering matrix and transmission vector from
///   the QuantumGraphObject and calculate the graph's eigenvalues, it
///   will probably be easier to create in one file, then print out the
///   quantum graph to a text file which downstream programs can use.
///
///  -Construction-
///   It currently makes most sense to declare an empty 
///   QuantumGraphObject and then fill it in with nodes and bonds. The
///   defualt constructor serves this purpose nicely:
///       QuantumGraphObject QGO;
///
///  -Adding Nodes-
///   In the base design here, a node must be created knowing exactly
///   how many bonds will be connected to it and the scattering matrix
///   that will apply to it. After one defines a scattering matrix sigma
///   as a gsl_matrix_complex* type, we can add a node to the graph 
///   object by calling
///       QGO.AddNode(sigma);
///   after which the node is stored in a vector of pointers as 
///   QGO.Node[0], if it was the first node added.
///
///  -Adding Bonds-
///   
///
///  -Updating Matrices-
///   When the QuantumGraphObject is complete, and you have verified
///   that each node is connected to a bond and each bond is connected
///   to two nodes in the way you desire, the internal matrices of the
///   QuantumGraph base class still need to change to reflect the new
///   structure, node scattering, and bond lengths. You can do this by
///   running the update function
///       QGO.UpdateQuantumGraph();
///   which will automatically create the scattering matrix and vector
///   of bond lengths for the entire quantum graph. Note that this is
///   a bond-to-bond scattering matrix, and the length vector is really
///   a representation of the transmission matrix for the graph's bonds.
///
///  -Saving the Graph
///   The QuantumGraphObject is not needed for future calculations, but
///   the data stored by the QuantumGraph object is. The QuantumGraph
///   class has an output method that will store the representation of
///   the graph in text format. Simply use
///       std::cout << QGO << std::endl;
///   and redirect the program output to a text file. This will store a
///   version of the graph that eigenvalue and eigenvector finding 
///   programs can use further along in the workflow.
///
class QuantumGraphObject: public QuantumGraph
{

  private:
    // Contains a list of pointers to nodes that make up the graph. The
    // QuantumGraphObject class is responsible for the memory for these
    // pointers. The memory is freed in the QuantumGraphObject 
    // destructor, and nowhere else.
    std::vector<QuantumGraphNode*> nodes_;


    // Contains a list of pointers to Undirected Bonds that make up the
    // graph. The QuantumGraphObject class is responsible for the memory
    // for these pointers. The memory is freed in the QuantumGraphObject
    // destructor, and nowhere else.
    //
    // Note that the QuantumGraph Object Class does not directly keep a
    // list of its constituent directed bonds, but only the undirected
    // bonds. Any interfacing with the directed bonds must be done
    // through the undirected bond class. This ensures that the two
    // directed bonds continue to be associated with each other when
    // something in the graph changes.
    std::vector<QuantumGraphUndirectedBond*> undirected_bonds_;


  public:
    // Defualt Constructor. Calls the QuantumGraph Constructor and 
    // sets the Nodes and Undirected Bond vectors to empty vectors.
    QuantumGraphObject();


    // Copy Constructor. Creates a copy of each node and undirected bond
    // in a seperate place in memory: a deep copy.
    QuantumGraphObject(const QuantumGraphObject&);


    // Destructor. Frees the memory for each Node and Undirected Bond.
    // Important, as these are dynamically assigned.
    ~QuantumGraphObject();



    // Creates a QuantumGraphNode Object with the argument as the node's
    // scattering matrix. This Node is dynamically created in memory and
    // the pointer is then stored in the vector of Node pointers.
    void AddNode(gsl_matrix_complex*);


    // 
    void Connect(unsigned int, unsigned int, gsl_complex);


    //
    void UpdateQuantumGraph();


    //
    int GetBondIndexFromPointer(QuantumGraphBond*);
    int GetNodeIndexFromPointer(QuantumGraphNode*);


};







////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
/// QuantumGraphNode Class
///
///   
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







////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// QuantumGraphBond Class
///
///   
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







////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// QuantumGraphUndirectedBond Class
///
///   
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






////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// END OF HEADER FILE
////////////////////////////////////////////////////////////////////////
#endif
