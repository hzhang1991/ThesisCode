#include "dog_math.h"
#include "dogdefs.h"
#include "mesh.h"

// -------------------------------------------------------------
// Routine for computing the L2-projection of an input function
// onto an orthonormal Legendre basis
// -------------------------------------------------------------
void L2Project_Unst(const int istart, 
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
            const dTensor2&,dTensor2&))
{

    // starting and ending indices 
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

    // number of quadrature points
    assert_ge(QuadOrder, 1);
    assert_le(QuadOrder, 5);

    // Number of quadrature points
    int mpoints;
    switch( QuadOrder )
    {
        case 1:
            mpoints = 1;
            break;

        case 2:
            mpoints = 3;
            break;

        case 3:
            mpoints = 6;
            break;

        case 4:
            mpoints = 12;
            break;

        case 5:	     
            mpoints = 16;
            break;
    }

    const int kmax = iMax(iMax(kmax_qin, kmax_auxin), kmax_fout);
    dTensor2  phi(mpoints, kmax); // Legendre basis (orthogonal)
    dTensor2 spts(mpoints, 2);    // List of quadrature points
    dTensor1 wgts(mpoints);       // List of quadrature weights

    void setQuadPoints_Unst(int QuadOrder, dTensor1& wgts, dTensor2& spts);
    setQuadPoints_Unst( QuadOrder, wgts, spts );

    // Evaluate the basis functions at each point
    void SetLegendreAtPoints_Unst(const dTensor2& spts, dTensor2& phi);
    SetLegendreAtPoints_Unst(spts, phi);

    // -------------------------------------------------------------
    // Loop over every grid cell indexed by user supplied parameters
    // described by istart...iend
    // -------------------------------------------------------------
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {      
        //find center of current cell
        const int    i1 = Mesh.get_tnode(i,1);
        const int    i2 = Mesh.get_tnode(i,2);
        const int    i3 = Mesh.get_tnode(i,3);

        // Corners:
        const double x1 = Mesh.get_node(i1,1);
        const double y1 = Mesh.get_node(i1,2);
        const double x2 = Mesh.get_node(i2,1);
        const double y2 = Mesh.get_node(i2,2);
        const double x3 = Mesh.get_node(i3,1);
        const double y3 = Mesh.get_node(i3,2);

        // Center of current cell:
        const double xc = (x1+x2+x3)/3.0;
        const double yc = (y1+y2+y3)/3.0;

        // Variables that need to be written to, and therefore are 
        // created for each thread
        dTensor2 xpts(mpoints,2);
        dTensor2 qvals(mpoints,meqn);
        dTensor2 auxvals(mpoints,maux);
        dTensor2 fvals(mpoints,mcomps_out);

        // Compute q, aux and fvals at each Gaussian Quadrature point
        // for this current cell indexed by (i,j)
        // Save results into dTensor2 qvals, auxvals and fvals.
        for (int m=1; m<=(mpoints); m++)
        {
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
        //
        // TODO - do we want to optimize this by looking into using transposes,
        // as has been done in the 2d/cart code? (5/14/2014) -DS
        for (int m1=1; m1<=mcomps_out; m1++)		
        for (int m2=1; m2<=kmax; m2++)
        {
            double tmp = 0.0;
            for (int k=1; k<=mpoints; k++)
            {
                tmp = tmp + wgts.get(k)*fvals.get(k,m1)*phi.get(k,m2);
            }
            fout->set(i,m1,m2, 2.0*tmp );
        }

    }

}
