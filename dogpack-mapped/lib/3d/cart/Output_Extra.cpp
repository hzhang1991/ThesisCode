#include "DogSolverCart3.h"
#include "DogStateCart3.h"
#include "DogParamsCart3.h"
#include "DogParams.h"
#include "tensors.h"

void Output_Extra(const DogSolverCart3& solver, int n)
{
  // Just in case you care to do something with q for your outputting
  // pleasure ... 
  const dTensorBC5& q = solver.get_state().get_q();
}
