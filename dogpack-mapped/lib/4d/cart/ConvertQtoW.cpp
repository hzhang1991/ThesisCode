#include "dogdefs.h"
#include "tensors.h"
#include "DogParams.h"

// Convert conserved variables to characteristic variables
void ConvertQtoW(int ixy, 
		 const dTensorBC5& aux, 
		 const dTensorBC5& qold, 
                 dTensorBC5& dwp, 
		 dTensorBC5& dwm, 
		 dTensorBC5& w_cent,
		 void (*ProjectLeftEig)(int,
					const dTensor1&,
					const dTensor1&,
					const dTensor2&,
					dTensor2&))
{
    printf("\n");
    printf("   ConvertQtoW has not yet been implemented in 4D \n");
    printf("\n");
    exit(1);
}
