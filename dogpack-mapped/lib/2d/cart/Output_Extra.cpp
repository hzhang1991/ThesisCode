#include <cmath>  // for fabs
#include "DogSolverCart2.h"
#include "DogStateCart2.h"
#include "DogParamsCart2.h"
#include "DogParams.h"
#include "tensors.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>


void Output_Extra(const DogSolverCart2& solver, int n)
{
    // Just in case you care to do something with q for your outputting
    // pleasure ... 
    const dTensorBC4& q = solver.get_state().get_q();
}
