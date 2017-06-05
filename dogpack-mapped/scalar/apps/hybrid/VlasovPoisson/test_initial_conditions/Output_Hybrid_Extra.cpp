#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "MonomialsToLegendre.h"
#include "mesh.h"
#include <sstream>
#include <string>
#include "assert.h"

/* 
 * This routine prints all the output stuff to file.  For this problem, we
 * wish to have a slice of the velocity space.  In this case, we cut data
 * along the line y=0, and save it as qslice[framenum].dat.
 * 
 */
void Output_Hybrid_Extra(const mesh& Mesh, const dTensorBC5& q, 
        double t, int nframe, string outputdir)
{ 

    const int       mx = q.getsize(1);   assert_eq(mx,dogParamsCart2.get_mx());
    const int       my = q.getsize(2);   assert_eq(my,dogParamsCart2.get_my());
    const int NumElems = q.getsize(3);
    const int     meqn = q.getsize(4);
    const int     kmax = q.getsize(5);
  
    int kmax2d;
    switch( dogParams.get_space_order() )
    {

        case 3:
            kmax2d = 6;
            break;

        case 2:
            kmax2d = 3;
            break;

        case 1:
            kmax2d = 1;
            break;
        default:
            unsupported_value_error( dogParams.get_space_order() );
    }

    // Load the slice data:

    // -- Load in the indices that we wish to save -- //
    // Alternatively, we could recompute these ... //
    dTensor1 BndList( NumElems );  // There are at most NumElems in the boundary ...
    string slice_index = "./matlab/slice_index.dat";
    FILE* file;
    if( (file = fopen(slice_index.c_str(), "r")) == NULL)
    { 
        printf("File slice_index.dat doesn't exists"); 
        printf("Please run FindSliceIndex() from the folder ./matlab/\n");
        return;
    }

    int bndy_index;
    int i = 0;
    while( fscanf(file, "%d", &bndy_index ) != EOF )
    {
        i += 1;
        BndList.set( i, bndy_index );
        printf("(i, bndy_index) = (%d, %d)\n", i, bndy_index );
    }
    const int NumTriOnBndy = i;
    fclose( file );

    // -- Open the file for writing to -- //
    ostringstream basename;
    basename << setfill('0') << setw(4) << nframe;
    string fnameslice  = outputdir+"/qslice"+basename.str()+".dat";
    file = fopen(fnameslice.c_str(),"w");

    // -- Print out all the fun information -- //
    fprintf(file,"%24.16e\n",t);
    for( int k=1; k<= kmax; k++             )
    for( int mt=1; mt <= NumTriOnBndy; mt++ )
    for( int j=1; j<=q.getsize(2); j++      )
    for( int i=1; i<=q.getsize(1); i++      )
    {
        double tmp = q.get(i, j, BndList.get(mt), 1, k );
        fprintf(file,"%24.16e\n",tmp);
    }
    fclose(file);

printf("Just printed stuff of size: kmax: %d NumTriOnBndy: %d my: %d mx: %d\n", kmax, NumTriOnBndy, my, mx );

}
