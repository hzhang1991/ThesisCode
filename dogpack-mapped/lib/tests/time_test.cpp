#include <assert.h>
#include <math.h>
#include "tensors.h"
#include "time_commands.h"

void time_test_test()
{
    // demonstration that dynamic memory allocation and
    // deallocation requires time on the order of acos (expensive)
    time_test(dTensor6 t(3,28,2,5,6,10));
    time_test(dTensor2 t(3,28));
    time_test(dTensor1 t(28));
    time_test(b=acos(a));
    time_test(b=sqrt(a));
    time_test(b=a*c);
    printf("\n");
}

int main()
{
    time_test_test();
}
