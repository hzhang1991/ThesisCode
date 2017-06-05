//#include <math.h>
#include <cmath>
#include <stdio.h>

int main(int argc, char** argv)
{

    double x = 2.404825557695773;
    printf("j0(%2.15e) = %2.15e\n", 0., j0(0.) );
    printf("j0(%2.15e) = %2.15e\n", x, j0(x) );

}
