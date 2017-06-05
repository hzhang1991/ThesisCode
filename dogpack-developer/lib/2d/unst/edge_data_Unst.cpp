#include <new>
#include "tensors.h"
#include "edge_data_Unst.h"

//using namespace std;  // not needed because we don't print anything here

// Constructor
edge_data_Unst::edge_data_Unst(const int NumEdgesIn)
{

    int i,m,k;

    NumEdges = NumEdgesIn;

    // Gauss-Legendre Quadrature
    phi_left  = new dTensor3(NumEdges,MAX_1D_PTS,MAX_2D_KMAX);
    phi_right = new dTensor3(NumEdges,MAX_1D_PTS,MAX_2D_KMAX);

    // 1D quadrature weights and points
    wgts1d = new dTensor1(MAX_1D_PTS);
    xpts1d = new dTensor1(MAX_1D_PTS);

    for (i=1; i<=NumEdges; i++)
    for (m=1; m<=MAX_1D_PTS; m++)
    for (k=1; k<=MAX_2D_KMAX; k++)
    {
        phi_left->set(i,m,k, 0.0 );
        phi_right->set(i,m,k, 0.0 );
    }

    for (m=1; m<=MAX_1D_PTS; m++)
    {
        wgts1d->set(m, 0.0 );
        xpts1d->set(m, 0.0 );
    }

    // Gauss-Lobatto Quadrature
    GL_phi_left  = new dTensor3(NumEdges,MAX_1D_PTS+1,MAX_2D_KMAX);
    GL_phi_right = new dTensor3(NumEdges,MAX_1D_PTS+1,MAX_2D_KMAX);

    GL_wgts1d = new dTensor1(MAX_1D_PTS+1);
    GL_xpts1d = new dTensor1(MAX_1D_PTS+1);

    for (i=1; i<=NumEdges; i++)
    for (m=1; m<=MAX_1D_PTS+1; m++)
    for (k=1; k<=MAX_2D_KMAX; k++)
    {
        GL_phi_left->set(i,m,k, 0.0 );
        GL_phi_right->set(i,m,k, 0.0 );
    }

    for (m=1; m<=MAX_1D_PTS+1; m++)
    {
        GL_wgts1d->set(m, 0.0 );
        GL_xpts1d->set(m, 0.0 );
    }

}

// Copy constructor
edge_data_Unst::edge_data_Unst(const edge_data_Unst& another_edge)
{

    int i,m,k;

    NumEdges = another_edge.NumEdges;

    for (i=1; i<=NumEdges; i++)
    for (m=1; m<=MAX_1D_PTS; m++)
    for (k=1; k<=MAX_2D_KMAX; k++)
    {
        phi_left->set(i,m,k, another_edge.phi_left->get(i,m,k) );
        phi_right->set(i,m,k, another_edge.phi_right->get(i,m,k) );
    }

    for (m=1; m<=MAX_1D_PTS; m++)
    {
        wgts1d->set(m, another_edge.wgts1d->get(m) );
        xpts1d->set(m, another_edge.xpts1d->get(m) );
    }

    for (i=1; i<=NumEdges; i++)
    for (m=1; m<=MAX_1D_PTS+1; m++)
    for (k=1; k<=MAX_2D_KMAX; k++)
    {
        GL_phi_left ->set(i,m,k, another_edge.GL_phi_left ->get(i,m,k) );
        GL_phi_right->set(i,m,k, another_edge.GL_phi_right->get(i,m,k) );
    }

    for (m=1; m<=MAX_1D_PTS+1; m++)
    {
        GL_wgts1d->set(m, another_edge.GL_wgts1d->get(m) );
        GL_xpts1d->set(m, another_edge.GL_xpts1d->get(m) );
    }
}

// Destructor
edge_data_Unst::~edge_data_Unst()
{

    // Gauss-Legendre Quadrature
    delete phi_left;
    delete phi_right;

    delete wgts1d;
    delete xpts1d;

    // Gauss-Lobatto Quadrature
    delete GL_phi_left;
    delete GL_phi_right;

    delete GL_wgts1d;
    delete GL_xpts1d;

}
