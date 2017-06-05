#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

// The purpose of this module is to implement various 1D quadrature rules that are
// required throughout DoGPack.  In most cases, building a multi-dimensional
// quadrature rule requires performing a tensor product of 1D quadrature
// weights and points.  This module describes the 1D rules that can be used to
// define higher dimensional quadrature rules.
//
// Two versions of each rule have been implemented.  If one is interested in
// performing numerical integration, then the version that uses points and
// weights should be called.  On the other hand, if one is only interested in
// constructing *points*, for example with applying a positivity-preserving
// limiter, then one should use the version that only sets the quadrature
// points.
//
// http://en.wikipedia.org/wiki/Gaussian_quadrature

// Gaussian quadrature (also known as Gauss-Legendre quadrature)
void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);
void setGaussPoints1d(dTensor1& x1d);

// Gauss-Lobatto quadrature
void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
void setGaussLobattoPoints1d( dTensor1& x1d);

#endif
