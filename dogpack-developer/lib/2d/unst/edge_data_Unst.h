#ifndef _EDGEDATA_UNST_H_
#define _EDGEDATA_UNST_H_

//#include<cmath>
//#include "tensors.h"
//using namespace std;
// --> moved to edge_data_Unst.cpp for faster compilation

class edge_data_Unst
{
    public:

        // Constructor
        edge_data_Unst(const int NumEdgesIn);

        // Copy constructor
        edge_data_Unst(const edge_data_Unst& another_edge);

        // Destructor
        ~edge_data_Unst();

        // Gauss-Legendre Quadrature
        dTensor3* phi_left;         //(NumEdges,5,15);
        dTensor3* phi_right;        //(NumEdges,5,15);

        dTensor1* wgts1d;           // length(5);
        dTensor1* xpts1d;           // length(5);

        // Gauss-Lobatto Quadrature
        dTensor3* GL_phi_left;      //(NumEdges,6,15);
        dTensor3* GL_phi_right;     //(NumEdges,6,15);

        dTensor1* GL_wgts1d;        // length(6);
        dTensor1* GL_xpts1d;        // length(6);

    private:

        int NumEdges;

        // Problem dimensions (these will need to be modified if we add
        // higher-order accuracy)
        static const int MAX_1D_PTS  = 5;
        static const int MAX_2D_KMAX = 15;

};

#endif
