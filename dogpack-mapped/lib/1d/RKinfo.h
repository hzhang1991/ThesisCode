#ifndef _RKINFO_H_
#define _RKINFO_H_

// Runge-Kutta information
struct RKinfo
{
  int mstage;
  int num_stages;

  dTensor1* alpha1;
  dTensor1* alpha2;
  dTensor1* beta;

  // These two are needed for 5th order Stepping
  dTensor2* gamma;
  dTensor1* delta;

};

#endif
