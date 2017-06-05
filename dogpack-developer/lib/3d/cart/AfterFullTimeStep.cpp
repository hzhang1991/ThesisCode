#include "dogdefs.h"
#include "DogSolverCart3.h"

// Function that is called after a full time step
// (i.e., after all stages are complete)
void AfterFullTimeStep(DogSolverCart3& solver)
{
  const double t = solver.get_state().get_time();
  dTensorBC5& aux = solver.fetch_state().fetch_aux();
  dTensorBC5& q = solver.fetch_state().fetch_q();

  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mz   = q.getsize(3);
  const int meqn = q.getsize(4);
  const int kmax = q.getsize(5);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(4);
}
