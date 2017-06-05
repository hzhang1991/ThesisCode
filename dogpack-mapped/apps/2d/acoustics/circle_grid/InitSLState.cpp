///////////////////////////////////////////////////////////////////////////////
// Dummy function here, applications who wish to override this are required to
// do so in their own application directory.
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "SLState.h"

// derivative of electric field, E, E_t
void ElectricField_repeat(const dTensor2& xpts, const dTensor2& qvals, 
        const dTensor2& auxvals, dTensor2& e, void* data)
{
    const int numpts=xpts.getsize(1);
    double x, v, t;

    // grab the current time and the next time
    t = *(double*)data;

    double E, Et;

    for (int i=1; i<=numpts; i++)
    {
        x = xpts.get(i,1);
        v = xpts.get(i,2);

        E    = 0.25 * (t -   t*t ) * cos( pi * x );
        Et   = 0.25 * (1.0-2.0*t ) * cos( pi * x );

        e.set(i, 1, E    );
        e.set(i, 2, Et   );

    }

}

void InitSLState( 
    const dTensorBC4& q, const dTensorBC4& aux, SL_state& sl_state )
{

    void ConvertQ1dToQ2d(const int &mopt,
            const dTensorBC3& qin, dTensorBC4& qout);
    void ConvertQ2dToQ1d(const int &mopt, int istart, int iend, 
                     const dTensorBC4& qin, dTensorBC3& qout);
    void IntegrateQ1dMoment1(const dTensorBC4& q2d, dTensorBC3& q1d);

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int maux = aux.getsize(2);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int mpoints = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
    const int kmax1d  = int(sqrt(mpoints));

    // cell widths are uniform.
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    const double t      = sl_state.t;
    const double dt     = sl_state.dt;
    const double tn     = sl_state.tn;

    // conversion from q2d(x) = q1d(x)
    dTensor1 k2d(5);
    k2d.set(5,14);
    k2d.set(4, 9);
    k2d.set(3, 5);
    k2d.set(2, 2);
    k2d.set(1, 1);

//  //////////////////////// Compute Electric Field E(t) ////////////////////
    dTensorBC3 Enow(mx, meqn, kmax1d, mbc,1);

    //////////// Compute Estar_t = rho_u - \int_{x,v} psi /////////////////////
    dTensorBC3 Enow_t(mx, meqn, kmax1d,mbc,1); 

    if ( dogParams.get_source_term()>0 )
    {
        dTensorBC4* ExactE;
        dTensorBC3* ExactE1d;

        double* t_ptr = new double;
        *t_ptr = tn;
        ExactE     = new dTensorBC4(mx,1,2,kmax,mbc);
        ExactE1d   = new dTensorBC3(mx,2,kmax1d,mbc);
	void L2Project_extra(const int istart, 
			     const int iend, 
			     const int jstart, 
			     const int jend,
			     const int QuadOrder, 
			     const int BasisOrder_qin,
			     const int BasisOrder_auxin,
			     const int BasisOrder_fout,		  
			     const dTensorBC4* qin,
			     const dTensorBC4* auxin, 
			     dTensorBC4* fout,
			     void (*Func)(const dTensor2&,const dTensor2&,
					  const dTensor2&,dTensor2&,void* data),
			     void* data);
	const int space_order = dogParams.get_space_order();
        L2Project_extra(1-mbc, mx+mbc, 1, 1, space_order, space_order, space_order, space_order,
			&q, &aux, ExactE, &ElectricField_repeat, (void*)t_ptr );
        ConvertQ2dToQ1d(1, 1-mbc, mx+mbc, *ExactE, *ExactE1d);

        int me = 1;
        for( int i=1-mbc; i <= mx+mbc; i++ )
        for( int k=1; k <=kmax1d; k++ )
        {
            // add in E_t as computed from source term //
            Enow.set(i,me,k, ExactE1d->get(i,1,k) );
            Enow_t.set(i,me,k, ExactE1d->get(i,2,k) );
        }
        delete ExactE;
        delete ExactE1d;
        delete t_ptr;
    }

    //////////// Combine E, Et, Estar and Estar_t and add into aux_extra //////

    // add in the weights to the appropriate spots
    for( int i = 1-mbc; i <= mx+mbc; i++ )
    for( int k = 1; k <= kmax1d; k ++ )
    {
        
        double En    = Enow.get(i,1,k);
        double En_t  = Enow_t.get(i,1,k);
        sl_state.aux1d->set(i,2,1, k, En );
        sl_state.aux1d->set(i,2,2, k, En_t );

   }

}// end of function InitSLState
///////////////////////////////////////////////////////////////////////////////
