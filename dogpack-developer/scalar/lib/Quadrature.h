#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

// This file is an identical duplicate of $DOGPACK/lib/Quadtrature.h
//
// (As of 4/16/2014) -DS

// Gaussian quadrature (also known as Gauss-Legendre quadrature)
void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);
void setGaussPoints1d(dTensor1& x1d);

// Gauss-Lobatto quadrature
void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
void setGaussLobattoPoints1d( dTensor1& x1d);

#endif
