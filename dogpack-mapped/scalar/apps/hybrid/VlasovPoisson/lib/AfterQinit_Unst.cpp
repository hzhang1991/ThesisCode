#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "mesh.h"
#include "VlasovParams.h"

#include "IniDocument.h"  // for initializing the vlasovParams

void AfterQinit_Unst(const mesh& Mesh, dTensorBC5& qnew)
{

    // initialize the vlasov parameters:
    vlasovParams.init();

    const int mx       = qnew.getsize(1);
    const int my       = qnew.getsize(2);
    const int NumElems = qnew.getsize(3);
    const int meqn     = qnew.getsize(4);
    const int kmax     = qnew.getsize(5);
    const int mbc      = qnew.getmbc();

    const int NumPhysElems = Mesh.get_NumPhysElems();
    // limit at all of the quadrature points
//  void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);
//  if( dogParams.using_moment_limiter() )
//  { ApplyPosLimiter(aux, qnew); }

    // Compute the integration constant (ion density); this will be used for all future calls
    // to ComputeElectricField
    vlasovParams.rho0 = 0.;
    vlasovParams.area = 0.;
    for ( int n=1; n<=NumPhysElems; n++ )
    {
        // TODO - is this the correct method for integrating over all physical
        // space?
        for ( int i=1; i<=mx; i++           )
        for ( int j=1; j<=my; j++           )
        { 
            vlasovParams.rho0 += 
                qnew.get(i,j,n,1,1)*Mesh.get_area_prim(n)*dogParamsCart2.get_prim_vol(); 
        }
        vlasovParams.area += Mesh.get_area_prim(n);
    }
    vlasovParams.rho0 = vlasovParams.rho0 / vlasovParams.area;

    printf("*** AfterQinit_Unst ***\n");
    printf("***    rho0 = %2.15e, area = %2.15e\n\n", vlasovParams.rho0, vlasovParams.area );
//  printf("diff = %2.15e\n", vlasovParams.area - pi*(1.0-0.6*0.6) );

    // Compute average momentum:
//  intConst = 0.0;
//  for(int i=1; i<=mx; i++)
//  { intConst += rho_u.get(i,1,1); }
//  intConst *= dogParamsCart2.get_dx()/(dogParamsCart2.get_xhigh()-dogParamsCart2.get_xlow());
//  vlasovParams.set_avg_momemtum(intConst);

//  printf("*** In AfterQinit, we are saving rho_u0 = %2.15e\n", vlasovParams.get_avg_momentum() );

}
