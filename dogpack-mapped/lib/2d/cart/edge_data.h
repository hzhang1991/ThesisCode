#ifndef _EDGE_DATA_H_
#define _EDGE_DATA_H_

#include <new>
#include <cmath>
#include "constants.h"
#include "tensors.h"
#include "assert.h"
#include "debug.h"
#include "DogParams.h"
#include "DogParamsCart2.h"

class dTensor1;
class dTensor2;
class dTensor3;

using namespace std;

class edge_data
{

    private:
        // disabled (could be implemented efficiently via tensor class copyfrom() method)
        edge_data(const edge_data& another_edge); // Copy constructor
        edge_data& operator=(const edge_data& in); // assignment operator
	int KMAX_MAX;

    public:
        edge_data(); // Constructor
        ~edge_data(); // Destructor
        void init(); // initialize

        // 1D Gaussian quadrature weights and points
        dTensor1* wgts1d; //(5);
        dTensor1* xpts1d; //(5);

        // Legendre basis functions
        dTensor2* phi_xl; //(5,15);
        dTensor2* phi_xr; //(5,15);
        dTensor2* phi_yl; //(5,15);
        dTensor2* phi_yr; //(5,15);

        // weights times Legendre basis functions
        dTensor2* wght_phi_xl; //(5,15);
        dTensor2* wght_phi_xr; //(5,15);
        dTensor2* wght_phi_yl; //(5,15);
        dTensor2* wght_phi_yr; //(5,15);
};

#endif
