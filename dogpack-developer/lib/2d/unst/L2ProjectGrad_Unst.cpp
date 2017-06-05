#include "dog_math.h"
#include "dogdefs.h"
#include "mesh.h"
#include "MonomialsToLegendre.h"

// -------------------------------------------------------------------------- //
// All-purpose routine for computing the L2-projection
// of various functions onto the gradient of the Legendre basis
//     (Unstructured grid version)
//
// See also: L2ProjectGrad.cpp (for the structured version)
// -------------------------------------------------------------------------- //
void L2ProjectGrad_Unst(const int istart, 
        const int iend, 
        const int QuadOrder, 
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,
        const mesh& Mesh, 
        const dTensor3* qin, 
        const dTensor3* auxin, 
        dTensor3* fout, 
        void (*Func)(const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor3&))
{

    // starting and ending indeces
    const int   NumElems = Mesh.get_NumElems();
    assert_ge(istart,1);
    assert_le(iend,NumElems);

    // qin variable
    assert_eq(NumElems,qin->getsize(1));
    const int     meqn = qin->getsize(2);
    const int kmax_qin = qin->getsize(3);
    assert_eq(kmax_qin,(BasisOrder_qin*(BasisOrder_qin+1))/2);

    // auxin variable
    assert_eq(NumElems,auxin->getsize(1));
    const int       maux = auxin->getsize(2);
    const int kmax_auxin = auxin->getsize(3);
    assert_eq(kmax_auxin,(BasisOrder_auxin*(BasisOrder_auxin+1))/2);

    // fout variables
    assert_eq(NumElems,fout->getsize(1));
    const int mcomps_out = fout->getsize(2);
    const int  kmax_fout = fout->getsize(3);
    assert_eq(kmax_fout,(BasisOrder_fout*(BasisOrder_fout+1))/2);

    // ---------------------------------------------------------------------- //
    // number of quadrature points
    // ---------------------------------------------------------------------- //
    assert_ge(QuadOrder,1);
    assert_le(QuadOrder,5);
    int mpoints;
    switch ( QuadOrder )
    {
        case 1:
            mpoints = 0;
            break;

        case 2:
            mpoints = 1;
            break;

        case 3:
            mpoints = 6;
            break;

        case 4:
            mpoints = 7;
            break;

        case 5:	     
            mpoints = 16;
            break;
    }

    // trivial case
    if ( QuadOrder==1 )
    {
#pragma omp parallel for
        for (int i=istart; i<=iend; i++)
        for (int m=1; m<=mcomps_out; m++) 
        for (int k=1; k<=kmax_fout; k++) 
        {  fout->set(i,m,k, 0.0 );  }
    }
    else
    {

        const int kmax = iMax(iMax(kmax_qin, kmax_auxin), kmax_fout);
        dTensor2    spts(mpoints, 2);
        dTensor1    wgts(mpoints);
        dTensor2      mu(mpoints, kmax);         // monomial basis (non-orthogonal)
        dTensor2     phi(mpoints, kmax);         // Legendre basis (orthogonal)
        dTensor2  phi_xi(mpoints, kmax_fout);    //  xi-derivative of Legendre basis (orthogonal)
        dTensor2 phi_eta(mpoints, kmax_fout);    // eta-derivative of Legendre basis (orthogonal)

        // Get the quadrature weights and points
        void setQuadPointsGrad_Unst(int QuadOrder, dTensor1& wgts, dTensor2& spts);
        setQuadPointsGrad_Unst(QuadOrder, wgts, spts);

        // Compute the Legendre polynomials evaluated at each point
        void SetLegendreAtPoints_Unst(const dTensor2& spts, dTensor2& phi);
        SetLegendreAtPoints_Unst(spts, phi);

        // Compute the partial derivatives of the canonical elements
        void SetLegendreGrad_Unst( const dTensor2& spts, dTensor2& phi_xi, dTensor2& phi_eta );
        SetLegendreGrad_Unst( spts, phi_xi, phi_eta );

        // ------------------------------------------------------------------ //
        // Loop over every grid cell indexed by user supplied parameters
        // described by istart...iend
        // ------------------------------------------------------------------ //
#pragma omp parallel for
        for (int i=istart; i<=iend; i++)
        {	  

            // local variables that are modified (set) per thread
            dTensor2    xpts(mpoints, 2);
            dTensor2   qvals(mpoints, meqn);
            dTensor2 auxvals(mpoints, maux);
            dTensor3   fvals(mpoints, mcomps_out, 2);

            // These need to be defined locally.  Each mesh element carries its
            // own change of basis matrix, so these need to be recomputed for
            // each element.  The canonical derivatives, phi_xi, and phi_eta can
            // be computed and shared for each element.
            dTensor2   phi_x(mpoints, kmax_fout);   //   x-derivative of Legendre basis (orthogonal)
            dTensor2   phi_y(mpoints, kmax_fout);   //   y-derivative of Legendre basis (orthogonal)

            // Find center of current cell
            const int i1    = Mesh.get_tnode(i,1);
            const int i2    = Mesh.get_tnode(i,2);
            const int i3    = Mesh.get_tnode(i,3);
            const double x1 = Mesh.get_node(i1,1);
            const double y1 = Mesh.get_node(i1,2);
            const double x2 = Mesh.get_node(i2,1);
            const double y2 = Mesh.get_node(i2,2);
            const double x3 = Mesh.get_node(i3,1);
            const double y3 = Mesh.get_node(i3,2);

            const double xc = (x1+x2+x3)/3.0;
            const double yc = (y1+y2+y3)/3.0;

            // Compute q, aux and fvals at each Gaussian Quadrature point
            // for this current cell indexed by (i,j)
            // Save results into dTensor2 qvals, auxvals and fvals.
            for (int m=1; m<=mpoints; m++)
            {
                // convert phi_xi and phi_eta derivatives
                // to phi_x and phi_y derivatives through Jacobian
                for (int k=1; k<=kmax_fout; k++)
                {
                    phi_x.set(m,k, Mesh.get_jmat(i,1,1)*phi_xi.get(m,k)
                                 + Mesh.get_jmat(i,1,2)*phi_eta.get(m,k) );
                    phi_y.set(m,k, Mesh.get_jmat(i,2,1)*phi_xi.get(m,k)
                                 + Mesh.get_jmat(i,2,2)*phi_eta.get(m,k) );
                }

                // point on the unit triangle
                const double s = spts.get(m,1);
                const double t = spts.get(m,2);

                // point on the physical triangle
                xpts.set(m,1, xc + (x2-x1)*s + (x3-x1)*t );
                xpts.set(m,2, yc + (y2-y1)*s + (y3-y1)*t );

                // Solution values (q) at each grid point
                for (int me=1; me<=meqn; me++)
                {
                    qvals.set(m,me, 0.0 );
                    for (int k=1; k<=kmax_qin; k++)
                    {
                        qvals.set(m,me, qvals.get(m,me) 
                                + phi.get(m,k) * qin->get(i,me,k) );
                    }
                }

                // Auxiliary values (aux) at each grid point
                for (int ma=1; ma<=maux; ma++)
                {
                    auxvals.set(m,ma, 0.0 );
                    for (int k=1; k<=kmax_auxin; k++)
                    {
                        auxvals.set(m,ma, auxvals.get(m,ma) 
                                + phi.get(m,k) * auxin->get(i,ma,k) );
                    }
                } 
            }

            // Call user-supplied function to set fvals
            Func(xpts,qvals,auxvals,fvals);

            // Evaluate integral on current cell (project onto Legendre basis) 
            // using Gaussian Quadrature for the integration
            for (int m1=1; m1<=mcomps_out; m1++)		
            for (int m2=1; m2<=kmax_fout; m2++)
            {
                double tmp = 0.0;
                for (int k=1; k<=mpoints; k++)
                {
                    tmp = tmp + wgts.get(k)*
                        ( fvals.get(k,m1,1)*phi_x.get(k,m2) +
                          fvals.get(k,m1,2)*phi_y.get(k,m2) );
                }
                fout->set(i,m1,m2, 2.0*tmp );
            }

        }
    }
}


// Similar function to L2ProjectGrad_Unst, but this time in place of using a user
// supplied function, we pass in Legendre weights for a function we already
// have the Legendre expansion of.
//
// Parameters:
// ----------
//
// The function we're projecting onto the basis functions is <F,G>:
//
//    F : F( 1:NumElems, 1:meqn, 1:kmax )
//    G : G( 1:NumElems, 1:meqn, 1:kmax )
//
// Returns:
// -------
//
// fout : the projected function of 
//
//        1/dA \iint \div( phi^{(k)} ) \cdot < F, G >.
//
//            
void L2ProjectGradAddLegendre_Unst(const int istart, const int iend, 
        const int QuadOrder, const mesh& Mesh, 
        const dTensor3* F, const dTensor3* G, dTensor3* fout )
{

    // starting and ending indeces
    const int   NumElems = Mesh.get_NumElems();
    assert_ge(istart,1);
    assert_le(iend,NumElems);

    // fout variables
    assert_eq(NumElems,fout->getsize(1));
    const int mcomps_out = fout->getsize(2);
    const int  kmax_fout = fout->getsize(3);

    // Number of equations
    const int meqn = F->getsize(2);   assert_eq( meqn, G->getsize(2) );

    // ---------------------------------------------------------------------- //
    // number of quadrature points
    // ---------------------------------------------------------------------- //
    assert_ge(QuadOrder,1);
    assert_le(QuadOrder,5);
    int mpoints;
    switch ( QuadOrder )
    {
        case 1:
            mpoints = 0;
            break;

        case 2:
            mpoints = 1;
            break;

        case 3:
            mpoints = 6;
            break;

        case 4:
            mpoints = 7;
            break;

        case 5:	     
            mpoints = 16;
            break;
    }

    // trivial case
    if ( QuadOrder==1 )
    {
#pragma omp parallel for
        for (int i = istart; i <= iend; i++)
        for (int m = 1; m <= mcomps_out; m++) 
        for (int k = 1; k <= kmax_fout; k++) 
        {  fout->set(i,m,k, 0.0 );  }
    }
    else
    {

        const int kmax = kmax_fout;
        dTensor2    spts(mpoints, 2);
        dTensor1    wgts(mpoints);
        dTensor2      mu(mpoints, kmax);         // monomial basis (non-orthogonal)
        dTensor2     phi(mpoints, kmax);         // Legendre basis (orthogonal)
        dTensor2  phi_xi(mpoints, kmax_fout);    //  xi-derivative of Legendre basis (orthogonal)
        dTensor2 phi_eta(mpoints, kmax_fout);    // eta-derivative of Legendre basis (orthogonal)

        // Get the quadrature weights and points
        void setQuadPointsGrad_Unst(int QuadOrder, dTensor1& wgts, dTensor2& spts);
        setQuadPointsGrad_Unst(QuadOrder, wgts, spts);

        // Compute the Legendre polynomials evaluated at each point
        void SetLegendreAtPoints_Unst(const dTensor2& spts, dTensor2& phi);
        SetLegendreAtPoints_Unst(spts, phi);

        // Compute the partial derivatives of the canonical elements
        void SetLegendreGrad_Unst( const dTensor2& spts, dTensor2& phi_xi, dTensor2& phi_eta );
        SetLegendreGrad_Unst( spts, phi_xi, phi_eta );

        // ------------------------------------------------------------------ //
        // Loop over every grid cell indexed by user supplied parameters
        // described by istart...iend
        // ------------------------------------------------------------------ //
#pragma omp parallel for
        for (int i=istart; i<=iend; i++)
        {	  

            // local variables that are modified (set) per thread
            dTensor2    xpts(mpoints, 2);
            dTensor3   fvals(mpoints, mcomps_out, 2);

            // These need to be defined locally.  Each mesh element carries its
            // own change of basis matrix, so these need to be recomputed for
            // each element.  The canonical derivatives, phi_xi, and phi_eta can
            // be computed and shared for each element.
            dTensor2   phi_x(mpoints, kmax_fout);   //   x-derivative of Legendre basis (orthogonal)
            dTensor2   phi_y(mpoints, kmax_fout);   //   y-derivative of Legendre basis (orthogonal)

            // Find center of current cell
            const int i1    = Mesh.get_tnode(i,1);
            const int i2    = Mesh.get_tnode(i,2);
            const int i3    = Mesh.get_tnode(i,3);
            const double x1 = Mesh.get_node(i1,1);
            const double y1 = Mesh.get_node(i1,2);
            const double x2 = Mesh.get_node(i2,1);
            const double y2 = Mesh.get_node(i2,2);
            const double x3 = Mesh.get_node(i3,1);
            const double y3 = Mesh.get_node(i3,2);

            const double xc = (x1+x2+x3)/3.0;
            const double yc = (y1+y2+y3)/3.0;

            // Compute q, aux and fvals at each Gaussian Quadrature point
            // for this current cell indexed by (i,j)
            // Save results into dTensor2 qvals, auxvals and fvals.
            for (int m=1; m<=mpoints; m++)
            {
                // convert phi_xi and phi_eta derivatives
                // to phi_x and phi_y derivatives through Jacobian
                for (int k=1; k<=kmax_fout; k++)
                {
                    phi_x.set(m,k, Mesh.get_jmat(i,1,1)*phi_xi.get(m,k)
                                 + Mesh.get_jmat(i,1,2)*phi_eta.get(m,k) );
                    phi_y.set(m,k, Mesh.get_jmat(i,2,1)*phi_xi.get(m,k)
                                 + Mesh.get_jmat(i,2,2)*phi_eta.get(m,k) );
                }

                // point on the unit triangle
                const double s = spts.get(m,1);
                const double t = spts.get(m,2);

                // point on the physical triangle
                xpts.set(m,1, xc + (x2-x1)*s + (x3-x1)*t );
                xpts.set(m,2, yc + (y2-y1)*s + (y3-y1)*t );

                // Solution values at each grid point
                for (int me=1; me <= meqn; me++)
                {
                    double val1 = F->get(i, me, 1);
                    double val2 = G->get(i, me, 1);
                    for (int k=2; k<=kmax; k++)
                    {
                        val1 += phi.get(m, k) * F->get(i, me, k);
                        val2 += phi.get(m, k) * G->get(i, me, k);
                    }
                    fvals.set(m, me, 1, val1 );
                    fvals.set(m, me, 2, val2 );
                }

            }

            // Call user-supplied function to set fvals
            //
            // NO NEED TO CALL THIS FUNCTION, WE ALREADY HAVE THE BASIS
            // FUNCTIONS!
            //
            // Func(xpts, qvals, auxvals, fvals);

            // Evaluate integral on current cell (project onto Legendre basis) 
            // using Gaussian Quadrature for the integration
            for (int me = 1; me <= mcomps_out; me++)		
            for (int k  = 1; k  <= kmax_fout;  k ++)
            {
                double tmp = 0.0;
                for (int m=1; m<=mpoints; m++)
                {
                    tmp += wgts.get(m) *
                        ( fvals.get(m, me, 1)*phi_x.get(m, k) +
                          fvals.get(m, me, 2)*phi_y.get(m, k) );
                }
                fout->set(i, me, k, fout->get(i, me, k) + 2.0*tmp );  // ADD to the RHS function
            }

        }
    }
    // printf("mcomps_out, kmax_fout = %d %d \n", mcomps_out, kmax_fout );

}
