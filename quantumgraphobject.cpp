
#include "quantumgraph.h"
#include "quantumgraphobject.h"

#include <iostream>

//
//
QuantumGraphObject::QuantumGraphObject(
  std::vector<QuantumGraphUndirectedBond*> UBondVec, 
  std::vector<QuantumGraphNode*> GraphNodeVec
)
{
  Nodes = GraphNodeVec;
  UndirectedBonds = UBondVec;
}


QuantumGraphObject::QuantumGraphObject()
{
  Nodes = {};
  UndirectedBonds = {};
}


void QuantumGraphObject::AddNode()
{
  Nodes.push_back(new QuantumGraphNode());
}


void QuantumGraphObject::AddNode(gsl_matrix_complex* S)
{
  Nodes.push_back(new QuantumGraphNode(S));
}


void QuantumGraphObject::Connect(QuantumGraphNode&, QuantumGraphNode&, 
                                 gsl_complex)
{
  
}


void QuantumGraphObject::Connect(unsigned int a, unsigned int b, 
                                 gsl_complex L)
{
  UndirectedBonds.push_back(new QuantumGraphUndirectedBond(Nodes[a], Nodes[b], L));
  Nodes[a]->ConnectToBond(UndirectedBonds.back());
  Nodes[b]->ConnectToBond(UndirectedBonds.back());
}


QuantumGraphObject::~QuantumGraphObject()
{
  std::clog << "Got here" << std::endl;
  for (unsigned int i=0; i<Nodes.size(); i++)
  {
    delete Nodes[i];
  }

  for (unsigned int i=0; i<UndirectedBonds.size(); i++)
  {
    delete UndirectedBonds[i];
  }
}


// 
//
//
void QuantumGraphObject::UpdateQuantumGraph()
{
  // For simplicity
  unsigned int numUBonds = UndirectedBonds.size();

  // Create somewhere to put the Scattering Matrix and the Length vector
  // When we're done.
  gsl_matrix_complex* newSMatrix 
      = gsl_matrix_complex_calloc(2 * numUBonds, 2 * numUBonds);
  gsl_vector_complex* newLVector
      = gsl_vector_complex_calloc(2 * numUBonds);


  // The 
  for (unsigned int bOut=0; bOut<numUBonds; bOut++)
  {
    std::vector<QuantumGraphUndirectedBond*> In;
    In = UndirectedBonds[bOut]->GetForwardBond()
                              ->GetStartNode()->GetConnectedBonds();
    for (unsigned int bIn=0; bIn<numUBonds; bIn++)
    {
      UndirectedBonds[bIn]->GetForwardBond();
    }

    In = UndirectedBonds[bOut]->GetBackwardBond()
                              ->GetStartNode()->GetConnectedBonds();
    for (unsigned int bIn=0; bIn<numUBonds; bIn++)
    {
      
    }
  }

  setGraph(newLVector, newSMatrix);

  gsl_matrix_complex_free(newSMatrix);
  gsl_vector_complex_free(newLVector);
}


//
//
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

std::vector<QuantumGraphUndirectedBond*> QuantumGraphNode::GetConnectedBonds()
{
  return ConnectedBonds;
}


//
//
QuantumGraphUndirectedBond::QuantumGraphUndirectedBond(
  QuantumGraphNode* ForwardStart,
  QuantumGraphNode* ForwardEnd, 
  gsl_complex LL
) :ForwardBond(ForwardStart, ForwardEnd, LL), 
   BackwardBond(ForwardEnd, ForwardStart, LL)
{
}


QuantumGraphBond* QuantumGraphUndirectedBond::GetForwardBond()
{
  return &ForwardBond;
}

QuantumGraphBond* QuantumGraphUndirectedBond::GetBackwardBond()
{
  return &BackwardBond;
}

//
//
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
