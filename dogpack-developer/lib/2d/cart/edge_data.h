#ifndef _EDGE_DATA_H_
#define _EDGE_DATA_H_

// Maximum size of polynomials considered
// For simplicity, we compute and save all the polynomials.  Because this class
// is instance once, and is the same for each element, this is not a problem.
#ifndef MAX_KMAX
#define MAX_KMAX 15
#endif

class dTensor1;
class dTensor2;
class dTensor3;
using namespace std;

class edge_data
{

    public:

        edge_data();        // Constructor
        ~edge_data();       // Destructor
        void init();        // initialize

        // 1D Gaussian quadrature weights and points.  Both of these have length
        // mpoints1d.
        dTensor1* wgts1d; 
        dTensor1* xpts1d;

        // Derivatives of the Legendre basis functions evaluated on each edge.  
        // Each of these will have size (1:mpoints1d, 1:MAX_KMAX ).
        dTensor2* phi_xl; 
        dTensor2* phi_xr; 
        dTensor2* phi_yl; 
        dTensor2* phi_yr;

        // Weights times derivatives of Legendre basis functions.  
        // Each of these will have size (1:mpoints1d, 1:MAX_KMAX ).
        //
        dTensor2* wght_phi_xl; 
        dTensor2* wght_phi_xr; 
        dTensor2* wght_phi_yl; 
        dTensor2* wght_phi_yr;

    private:

        // disabled (could be implemented efficiently via tensor class copyfrom() method)
        //
        // What do these operations do? -DS
        edge_data(const edge_data& another_edge);       // Copy constructor
        edge_data& operator=(const edge_data& in);      // assignment operator

        // Maximum polynomial order considered.  When creating an instance of
        // this class, we actually save enough information to go up to
        // 5th-order.
        // MOVED TO A MACRO, BECAUSE THIS GETS HARD CODED TO 15 ANYWAY (-DS).
//      int MAX_KMAX;  // set to 15 in constructor

};

#endif
