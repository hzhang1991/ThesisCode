#ifndef EDGEDATA_UNST_H
#define EDGEDATA_UNST_H

#include<cmath>
#include "tensors.h"
using namespace std;

class edge_data_Unst
{
 public:
  edge_data_Unst(const int NumEdgesIn);
  // Constructor
    
  edge_data_Unst(const edge_data_Unst& another_edge);
  // Copy constructor
    
  ~edge_data_Unst();
    // Destructor
  
  // Gauss-Legendre Quadrature
  dTensor3* phi_left; //(NumEdges,5,15);
  dTensor3* phi_right; //(NumEdges,5,15);
  
  dTensor1* wgts1d; //(5);
  dTensor1* xpts1d; //(5);

  // Gauss-Lobatto Quadrature
  dTensor3* GL_phi_left; //(NumEdges,6,15);
  dTensor3* GL_phi_right; //(NumEdges,6,15);
  
  dTensor1* GL_wgts1d; //(6);
  dTensor1* GL_xpts1d; //(6);
  
 private:
  int NumEdges;
};

#endif
