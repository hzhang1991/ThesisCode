#include "DogSolverCart4.h"
#include "DogStateCart4.h"
#include "DogParamsCart4.h"
#include "DogParams.h"
#include "tensors.h"

void Output_Extra(const DogSolverCart4& solver, int n)
{
  // Just in case you care to do something with q for your outputting
  // pleasure ... 
  const dTensorBC6& q = solver.get_state().get_q();
}
