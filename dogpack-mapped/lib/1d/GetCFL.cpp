#include <iostream>
#include "dog_math.h"
#include "tensors.h"
#include "stdlib.h"
using namespace std;

double GetCFL(double dt,double dtmax,
        const dTensor1& prim_vol,
        const int method[],
        const dTensorBC3& aux,
        const dTensorBC1& smax)
{
    const int melems =aux.getsize(1);
    double cfl=-100.0;

    if (dt>dtmax)
    {
        cout << endl;
        cout << " Error: dt is out of bounds ... " << endl;
        cout << "     dt = " << dt << endl;
        cout << endl;
        exit(1);
    }

// This can't be done in parallel!!!  (-DS)
//#pragma omp parallel for
    for (int j=1; j<=melems; j++)
    {
        cfl = Max(dt*smax.get(j)/prim_vol.get(j),cfl);  
    }

    if (cfl>1.0e8)
    {
        cout << endl;
        cout << " Error: CFL number is out of bounds ... " << endl;
        cout << "     CFL = " << cfl << endl;
        cout << endl;
        exit(1);
    }

    return cfl;
}
