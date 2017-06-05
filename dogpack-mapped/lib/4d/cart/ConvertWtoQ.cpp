#include "dogdefs.h"
#include "dog_math.h"

// Convert characteristic variables to conserved variables
void ConvertWtoQ(int ixy, 
		 const dTensorBC5& aux, 
		 const dTensorBC5& qold,
                 const dTensorBC5& win, 
		 dTensorBC5& qout,
		 void (*ProjectRightEig)(int,
					 const dTensor1&,
					 const dTensor1&,
					 const dTensor2&,
					 dTensor2&))
{
    printf("\n");
    printf("   ConvertWtoQ has not yet been implemented in 4D \n");
    printf("\n");
    exit(1);
}
