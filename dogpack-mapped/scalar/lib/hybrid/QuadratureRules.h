#ifndef _QUADRATURERULES_H_
#define _QUADRATURERULES_H_

#include "dogdefs.h"
#include <cmath>    
#include "MonomialsToLegendre.h"

class QuadratureRules
{

    private:

        // Quadrature points for cartesian, and unstructured grids:
        //
        // Their sizes are given by:
        //    spts( 1:num_quad_points, 1:2 )
        //
        dTensor2* spts_cart;
        dTensor2* spts_unst;

        // Quadrature weights for cartesian, and unstructured grids:
        //
        // Their sizes are given by:
        //    wgt( 1:num_quad_points )
        //  
        dTensor1* wgt_cart;
        dTensor1* wgt_unst;

        // Basis functions for cartestian, and unstructured grids:
        //
        // In both cases, their sizes are given by:
        //    phi( 1:num_quad_points, 1:kmax2d ).
        //
        // Note: there are a different number of quadrature points for
        // 2dunst as compared with 2dcart.
        dTensor2* phi_cart;
        dTensor2* phi_unst;

        // quadrature weights, points, and basis functions evaluated at this
        // list.  The format here allows the user to access the unstructured,
        // as well as the structured parts of the quadrature rules.
        //
        // Sizes are given by:
        //
        //  spts( 1:mpoints_unst, 1:mpoints_cart, 1:4    )
        //   wgt( 1:mpoints_unst, 1:mpoints_cart )
        //   phi( 1:mpoints_unst, 1:mpoints_cart, 1:kmax )
        //
        dTensor3* spts;
        dTensor2* wgt;
        dTensor3* phi;

        // sorder - spatial order of accuracy.
        //
        // kmax2d - number of 2D polynomial basis functions used.  Each
        // structured, and unstructured routines use the exact same number
        // here.
        //
        // kmax - number of basis functions used for the 4D representation.
        //
        int sorder, kmax2d, kmax;

    public:

        // The initial values for all private variables are set to zero.
        void init( int sorder );

        // Constructor:
        QuadratureRules(): wgt_unst (NULL), wgt_cart (NULL), wgt(NULL),
                           spts_unst(NULL), spts_cart(NULL), spts(NULL),
                           phi_unst (NULL), phi_cart (NULL), phi(NULL)
                           { sorder = -1; kmax2d = 0; kmax = 0; }

        // Desctructor:
        ~QuadratureRules();

        // Methods for evaluating the basis functions at each quadrature point
        //
        // G_unst( 1:meqn, 1:kmax2d )
        //
        // G_cart( 1:meqn, 1:kmax2d )
        //
        void EvaluateBasisUnst( const dTensor2& G_unst,  dTensor2& gvals );
        void EvaluateBasisCart( const dTensor2& G_cart,  dTensor2& gvals );

        // Accessors (declared in the header to avoid function call overhead)
        //
        const dTensor2& get_phi_unst()  const{ return *phi_unst;  } 
        const dTensor2& get_phi_cart()  const{ return *phi_cart;  } 

        const dTensor2& get_spts_unst() const{ return *spts_unst; } 
        const dTensor2& get_spts_cart() const{ return *spts_cart; } 

        const dTensor1& get_wgt_cart()  const{ return *wgt_cart; }
        const dTensor1& get_wgt_unst()  const{ return *wgt_unst; }

        const dTensor3& get_phi()       const{ return *phi;  } 
        const dTensor3& get_spts()      const{ return *spts; } 
        const dTensor2& get_wgt()       const{ return *wgt; }

        const int& get_kmax()           const{ return kmax; }
        const int& get_kmax2d()         const{ return kmax2d; }

};

#endif
