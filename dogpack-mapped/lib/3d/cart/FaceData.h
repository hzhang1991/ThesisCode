#ifndef _FACE_DATA_H_
#define _FACE_DATA_H_

#include <new>
#include <cmath>
#include "constants.h"
#include "tensors.h"
#include "assert.h"
#include "debug.h"
#include "DogParams.h"
#include "DogParamsCart3.h"

class dTensor1;
class dTensor2;
class dTensor3;

using namespace std;

class FaceData
{
 private:
  // disabled (could be implemented efficiently via tensor class copyfrom() method)
  FaceData(const FaceData& another_face); // Copy constructor
  FaceData& operator=(const FaceData& in); // assignment operator
  int KMAX_MAX;
  
 public:
  FaceData(); // Constructor
  ~FaceData(); // Destructor
  void init(); // initialize
  
  // 2D Gaussian quadrature weights and points
  dTensor1* wgts2d; //(25);
  dTensor2* xpts2d; //(25,2);
  
  // Legendre basis functions
  dTensor2* phi_xl; //(25,20);
  dTensor2* phi_xr; //(25,20);
  dTensor2* phi_yl; //(25,20);
  dTensor2* phi_yr; //(25,20);
  dTensor2* phi_zl; //(25,20);
  dTensor2* phi_zr; //(25,20);
  
  // weights times Legendre basis functions
  dTensor2* wght_phi_xl; //(25,20);
  dTensor2* wght_phi_xr; //(25,20);
  dTensor2* wght_phi_yl; //(25,20);
  dTensor2* wght_phi_yr; //(25,20);
  dTensor2* wght_phi_zl; //(25,20);
  dTensor2* wght_phi_zr; //(25,20);
};

#endif
