// stuff copied from apps/2d/VlasovPoisson1d/lib/Output.cpp
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "DogParams.h"
#include "DogParamsCart2.h"
#include "dogdefs.h"
#include "mesh.h"
#include "constants.h"

// Save the time, and print the L2-norm of the electric field
//
// TODO - is this the spot where we'd like to print moments?  Do we need
// moments after each time step or for each frame?
//
void PrintElectricField( const double t, const mesh& Mesh, const dTensorBC5& q, 
        const dTensor2& E1, const dTensor2& E2 )
{

    const int mx        = q.getsize(1);
    const int my        = q.getsize(2);
    const int NumElems  = q.getsize(3);
    const int meqn      = q.getsize(4);
    const int kmax      = q.getsize(5);
    const int mbc       = q.getmbc();

    const int NumPhysElems = Mesh.get_NumPhysElems();

    string fname1 = string(dogParams.get_outputdir())+"/conservation.dat";
    string fname2 = string(dogParams.get_outputdir())+"/Efield.dat";

    ofstream write_file1,write_file2;

    if( fabs(t) < EPSILON )
    {
        write_file1.open(fname1.c_str(), ofstream::out);
        write_file2.open(fname2.c_str(), ofstream::out);
    }
    else
    {
        write_file1.open(fname1.c_str(), ofstream::app);
        write_file2.open(fname2.c_str(), ofstream::app);
    }

    // -----------------
    // CONSERVATION
    // -----------------
    dTensor1 qsum(meqn);
    if (dogParams.get_mcapa()<1) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m,0.0);

            for (int i=1; i<=mx; i++)	    
            for (int j=1; j<=my; j++)
            for (int n=1; n<=NumPhysElems; n++)
            {
                qsum.set(m, qsum.get(m) + Mesh.get_area_prim(n)*dogParamsCart2.get_prim_vol()*q.get(i,j,n,m,1) );
            }

        }
    }

    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        if (fabs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
        write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
    write_file1 << endl;
    write_file1.close();

    ///////////////////////////////////////////////////////////////////////////
    // Conserved Vlasov-Poisson Quantities                                   //
    //        ||f||_1, ||f||_2, Energy, Entropy                              //
    ///////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // Electric Field 
    const int mcons = 4;

// TODO - do we want to compute the extra moments here or later?
////dTensorBC4 ConsData2d(mx, my, mcons, kmax,mbc); 
////const int space_order = dogParams.get_space_order();
////L2Project(1, mx, 1, my, space_order, space_order, space_order,
////      space_order, &q, &aux, &ConsData2d, &ConservedFunc );

    // electric field, entropy and energy
    dTensor1 ConsData(mcons); 
    ConsData.setall(0.);

    // Compute the integral of all the 2d quantities:
//  for (int m=1; m<=mcons; m++)
//  {
//      for (int i=1; i<=mx; i++)	    
//      for (int j=1; j<=my; j++)
//      {
//          ConsData.set(m, ConsData.get(m) + dx*dy*ConsData2d.get(i,j,m,1) );
//      }
//  }

    // Compute the integral of all the 1d quanties:
    //        \int E^2 \ dx
    double Esqd = 0.0;
    for(int n=1; n <= NumElems; n++ )
    for(int k=1; k <= E1.getsize(2); k++ )
    { 
        Esqd += Mesh.get_area_prim(n) * ( 
                pow( E1.get(n,k), 2 ) + pow(E2.get(n,k), 2 ) 
            ); 
    }

//printf("Esqd = %f\n", Esqd );
//exit(1);

    // Add in E to the total Energy:
    // ConsData.set(3, ConsData.get(3) + 0.5 * E2 );
    ConsData.set(3, 0. );

    // print time:
    write_file2 << setprecision(16);
    write_file2 << setw(24) << scientific << t << " ";

    // print ||E||_2 first:
    write_file2 << setw(24) << scientific << sqrt( Esqd ) << " ";

    // print all the other conserved quantities:
    for (int m=1; m<=mcons; m++)
    {
        write_file2 << setw(24) << scientific << ConsData.get(m) << " ";
    }
    write_file2 << endl;

    write_file2.close();

}
