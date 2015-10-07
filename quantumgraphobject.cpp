
#include "quantumgraph.h"
#include "quantumgraphobject.h"

#include <iostream>
#include <vector>

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
  for (unsigned int i=0; i<Nodes.size(); i++)
  {
    delete Nodes[i];
  }

  for (unsigned int i=0; i<UndirectedBonds.size(); i++)
  {
    delete UndirectedBonds[i];
  }
}


int QuantumGraphObject::GetBondIndexFromPointer(QuantumGraphBond* QGBP)
{
  for (unsigned int ub=0; ub<UndirectedBonds.size(); ub++)
  {
    if ((UndirectedBonds[ub]->GetForwardBond()) == QGBP)
    {
      return 2 * ub;
    }
    else if ((UndirectedBonds[ub]->GetBackwardBond()) == QGBP)
    {
      return 2 * ub + 1;
    }
  }

  return -1;
}


int QuantumGraphObject::GetNodeIndexFromPointer(QuantumGraphNode* QGNP)
{
  for (unsigned int n=0; n<Nodes.size(); n++)
  {
    if (Nodes[n] == QGNP)
    {
      return n;
    }
  }

  return -1;
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
  for (unsigned int b=0; b<numUBonds; b++)
  {
    //debug
    //std::clog << std::endl;
   // std::clog << "b = " << b << std::endl;
   // std::clog << "   ForwardBond" << std::endl;
    //Calculate for the Forward Bond
    int i = 2*b;
    gsl_vector_complex_set(newLVector, i, UndirectedBonds[b]->GetForwardBond()->GetBondLength());

    int l;
    QuantumGraphBond* ForwardBond = UndirectedBonds[b]->GetForwardBond();
    QuantumGraphNode* FStartNodeP = ForwardBond->GetStartNode();
    std::vector<QuantumGraphBond*> FIn = FStartNodeP->GetIncomingBonds();

    l = FStartNodeP->GetBondIndexAtNodeFromPointer(ForwardBond);

      // std::clog << "     Comes from node "<< GetNodeIndexFromPointer(FStartNodeP) << std::endl;
    for (unsigned int m=0; m<FIn.size(); m++)
    {
       int j = GetBondIndexFromPointer(FIn[m]);
       //std::clog << "     S Element =  ("<< i << "," << j << ")" << std::endl;

      // std::clog << "     Node S  Element =  (" << l << "," << m << ")" << std::endl;

       gsl_complex Sij = FStartNodeP->GetMatrixElement(l, m);
       gsl_matrix_complex_set(newSMatrix, i, j, Sij);
    }


    //debug
    //std::clog << "   BackwardBond" << std::endl;
    // Now For the Backward Bond
    i = (2*b)+1;

    gsl_vector_complex_set(newLVector, i, UndirectedBonds[b]->GetBackwardBond()->GetBondLength());

    QuantumGraphBond* BackwardBond = UndirectedBonds[b]->GetBackwardBond();
    QuantumGraphNode* BStartNodeP = BackwardBond->GetStartNode();
    std::vector<QuantumGraphBond*> BIn = BStartNodeP->GetIncomingBonds();

    l = BStartNodeP->GetBondIndexAtNodeFromPointer(BackwardBond);
   //    std::clog << "     Comes from node "<< GetNodeIndexFromPointer(BStartNodeP) << std::endl;
    for (unsigned int m=0; m<BIn.size(); m++)
    {
       int j = GetBondIndexFromPointer(BIn[m]);
      // std::clog << "     S Element =  ("<< i << "," << j << ")" << std::endl;

       //std::clog << "     Node S  Element =  (" << l << "," << m << ")" << std::endl;

       gsl_complex Sij = BStartNodeP->GetMatrixElement(l, m);
       gsl_matrix_complex_set(newSMatrix, i, j, Sij);
    }
    //*/

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


int QuantumGraphNode::GetBondIndexAtNodeFromPointer(QuantumGraphBond* QGBP)
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

std::vector<QuantumGraphUndirectedBond*> QuantumGraphNode::GetConnectedUBonds()
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

QuantumGraphNode* QuantumGraphBond::GetEndNode()
{
  return EndNode;
}

gsl_complex QuantumGraphBond::GetBondLength()
{
  return ComplexLength;
}


