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

    // Compute the integration constant (ion density); this will be used for all future calls
    // to ComputeElectricField
    vlasovParams.rho0 = 0.;
    vlasovParams.area = 0.;
    for ( int n=1; n<=NumPhysElems; n++ )
    {

        for( int i=1; i<=mx; i++ )
        for( int j=1; j<=my; j++ )
        { 
            vlasovParams.rho0 += 
                qnew.get(i,j,n,1,1)*Mesh.get_area_prim(n)*dogParamsCart2.get_prim_vol(); 
        }
        vlasovParams.area += Mesh.get_area_prim(n);
    }
    vlasovParams.rho0 = vlasovParams.rho0 / vlasovParams.area;

    printf("Mesh.get_area_prim(1) = %f; dogParamsCart2.get_prim_vol() = %f\n", Mesh.get_area_prim(1), dogParamsCart2.get_prim_vol() );
    printf("*** AfterQinit_Unst ***\n");
    printf("***    rho0 = %2.15e, area = %2.15e\n\n", vlasovParams.rho0, vlasovParams.area );

    // For testing Landau Damping problem
//  printf("       diff in area    = %2.15e\n", vlasovParams.area - (4.0*pi)*(4.0*pi) );
//  printf("       diff in density = %2.15e\n", vlasovParams.rho0 - 1.0 );

}
