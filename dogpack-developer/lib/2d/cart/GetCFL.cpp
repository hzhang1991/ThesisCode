#include "tensors.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogSolverCart2.h"

double DogSolverCart2::GetCFL(double dt) const
{
  const dTensorBC3& smax = get_smax();
  const int mx = smax.getsize(1);
  const int my = smax.getsize(2);
  
  const double prim_vol = dogParamsCart2.get_prim_vol();
  const double dt_over_prim_vol = dt/prim_vol;
  dTensor1 cflx(mx); // parallelization requires memory allocation
  cflx.setall(0.);
  
#pragma omp parallel for
  for (int i=1; i<=mx; i++)
    {
      double cfl_i=0.;
      for (int j=1; j<=my; j++)
	{
          const double tmp = Max(smax.get(i,j,1),smax.get(i,j,2));
          const double local_cfl = tmp*dt_over_prim_vol;
          cfl_i = Max(local_cfl, cfl_i);
	}
      cflx.set(i, cfl_i);
    }
  
  double cfl=0.;
  for (int i=1; i<=mx; i++)
    {
      cfl = Max(cflx.get(i), cfl);
    }
  
  return cfl;
}
