#include "../defs.h"

// Function that is called before each time step
void BeforeStep_sub(dTensor2 node, dTensorBC3& aux, dTensorBC3& q)
{
    int j;
    double x; 
    int melems = q.getsize(1);
    int   meqn = q.getsize(2);
    int   kmax = q.getsize(3);
    int   maux = aux.getsize(2);
    
}
