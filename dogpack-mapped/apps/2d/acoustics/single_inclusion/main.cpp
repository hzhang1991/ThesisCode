#include "AppSolver.h"

// =========================================================================
//
//  Copyright J.A. Rossmanith
//
//  This software is made available for research and instructional use only.
//  You may copy and use this software without charge for these non-commercial
//  purposes, provided that the copyright notice and associated text is
//  reproduced on all copies.  For all other uses (including distribution of
//  modified versions), please contact the author at the address given below.
//
//  *** This software is made available "as is" without any assurance that it
//  *** will work for your purposes.  The software may in fact have defects, so
//  *** use the software at your own risk.
//
//  -------------------------------
//  DoGPack
//  -------------------------------
//
//    Lead Developer:  
//             James Rossmanith
//             Iowa State University
//             Department of Mathematics
//             396 Carver Hall
//             Ames, IA 50010
//             rossmani@iastate.edu
// =========================================================================
//
// NOTE: To define your solver, you can modify the following files:
// 1. parameters.ini -- basic data file, can modify
//       # of grid points, time step, order of
//       accuracy in both space and time, etc...
// 2. AppSolver.h and AppSolver.cpp,
// 3. SetBndValues.cpp -- boundary conditions file, and
// 4. The following linking callback files:
//    1. called via L2Project:
//       1. QinitFunc.cpp -- initial condition file
//       2. AuxFunc.cpp -- auxiliary variable initial condition file
//       3. SourceTermFunc.cpp -- source term file
//    2. FluxFunc.cpp -- flux function file
//       (called via L2ProjectGrad and from RiemannSolve),
//    3. SetWaveSpd.cpp -- eigenvalues of flux Jacobian file
//       (called from RiemannSolve)
//    4. called from ApplyLimiter:
//       1. ProjectLeftEig.cpp -- left eigenvectors of flux Jacobian
//       2. ProjectLeftEig.cpp -- right eigenvectors of flux Jacobian
//
int main(int argc, char* argv[])
{
  AppSolver appSolver;
  return appSolver.main(argc, argv);
}
