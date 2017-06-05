#ifndef DIFFUSION_H
#define DIFFUSION_H
#include "constants.h"

const double phiv[5] = {1.0, sq3, sq5, sq7, 3.0};
const double Smat[5][5] = {{0.0, 0.0, 0.0, 0.0, 0.0},
                          {2.0*sq3, 0.0, 0.0, 0.0, 0.0},
                          {0.0, 2.0*sq3*sq5, 0.0, 0.0, 0.0},
                          {2.0*sq7, 0.0, 2.0*sq5*sq7, 0.0, 0.0},
                          {0.0, 6.0*sq3, 0.0, 6.0*sq7, 0.0}};
const double signv[5] = {-1.0, 1.0, -1.0, 1.0, -1.0};

#endif
