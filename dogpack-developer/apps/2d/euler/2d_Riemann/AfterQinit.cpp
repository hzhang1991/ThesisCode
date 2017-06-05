/*
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "DogParams.h"
*/

#include "DogSolverCart2.h"
#include <string.h>          // For strcpy and strcat (some compilers don't need this)
#include "EulerParams.h"

// Function that is called after initial condition
void AfterQinit(DogSolverCart2& solver)
{

//  DogStateCart2& state = solver.fetch_state();
//  dTensorBC4& aux = state.fetch_aux();
//  dTensorBC4& q   = state.fetch_q();
//  const int mx    = q.getsize(1);
//  const int my    = q.getsize(2);
//  const int meqn  = q.getsize(3);
//  const int kmax  = q.getsize(4);
//  const int mbc   = q.getmbc();
//  const int maux  = aux.getsize(3);  
//  const int space_order = dogParams.get_space_order();

    // Output parameters to file in outputdir
    char eulerhelp[200];
    strcpy( eulerhelp, solver.get_outputdir() );
    strcat( eulerhelp, "/eulerhelp.dat");
    eulerParams.write_eulerhelp( eulerhelp );

}
