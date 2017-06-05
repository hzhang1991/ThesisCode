#include "edge_data_Unst.h"
#include<cmath>
#include<new>
#include "tensors.h"
using namespace std;

edge_data_Unst::edge_data_Unst(const int NumEdgesIn)
// Constructor
{
  int i,m,k;

  NumEdges = NumEdgesIn;

  // Gauss-Legendre Quadrature
  phi_left = new dTensor3(NumEdges,5,15);
  phi_right = new dTensor3(NumEdges,5,15);
  
  wgts1d = new dTensor1(5);
  xpts1d = new dTensor1(5);
  
  for (i=1; i<=NumEdges; i++)
    for (m=1; m<=5; m++)
      for (k=1; k<=15; k++)
	{
	  phi_left->set(i,m,k, 0.0 );
	  phi_right->set(i,m,k, 0.0 );
	}

  for (m=1; m<=5; m++)
    {
      wgts1d->set(m, 0.0 );
      xpts1d->set(m, 0.0 );
    }

  // Gauss-Lobatto Quadrature
  GL_phi_left  = new dTensor3(NumEdges,6,15);
  GL_phi_right = new dTensor3(NumEdges,6,15);

  GL_wgts1d = new dTensor1(6);
  GL_xpts1d = new dTensor1(6);
  
  for (i=1; i<=NumEdges; i++)
    for (m=1; m<=6; m++)
      for (k=1; k<=15; k++)
	{
	  GL_phi_left->set(i,m,k, 0.0 );
	  GL_phi_right->set(i,m,k, 0.0 );
	}

  for (m=1; m<=6; m++)
    {
      GL_wgts1d->set(m, 0.0 );
      GL_xpts1d->set(m, 0.0 );
    }

}

edge_data_Unst::edge_data_Unst(const edge_data_Unst& another_edge)
// Copy constructor
{
  int i,m,k;

  NumEdges = another_edge.NumEdges;

  for (i=1; i<=NumEdges; i++)
    for (m=1; m<=5; m++)
      for (k=1; k<=15; k++)
	{
	  phi_left->set(i,m,k, another_edge.phi_left->get(i,m,k) );
	  phi_right->set(i,m,k, another_edge.phi_right->get(i,m,k) );
	}

  for (m=1; m<=5; m++)
    {
      wgts1d->set(m, another_edge.wgts1d->get(m) );
      xpts1d->set(m, another_edge.xpts1d->get(m) );
    }

  for (i=1; i<=NumEdges; i++)
    for (m=1; m<=6; m++)
      for (k=1; k<=15; k++)
	{
	  GL_phi_left->set(i,m,k, another_edge.GL_phi_left->get(i,m,k) );
	  GL_phi_right->set(i,m,k, another_edge.GL_phi_right->get(i,m,k) );
	}

  for (m=1; m<=6; m++)
    {
      GL_wgts1d->set(m, another_edge.GL_wgts1d->get(m) );
      GL_xpts1d->set(m, another_edge.GL_xpts1d->get(m) );
    }
}

edge_data_Unst::~edge_data_Unst()
// Destructor
{
}
