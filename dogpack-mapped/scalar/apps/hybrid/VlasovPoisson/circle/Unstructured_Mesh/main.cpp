#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "meshdefs.h"

// =========================================================================
//
//  Copyright M. Elsey and J.A. Rossmanith
//
//  This software is made available for research and instructional use only.
//  You may copy and use this software without charge for these non-commercial
//  purposes, provided that the copyright notice and associated text is
//  reproduced on all copies.  For all other uses (including distribution of
//  modified versions), please contact the author at the address given below.
//
//  *** This software is made available "as is" without any assurance that it
//  *** will work for your purposes.  The software may in fact have defects, so
//  *** use the software at your own risk.
//
//    Authors:  Matthew Elsey
//              Department of Mathematics
//              University of Michigan
//              Ann Arbor, MI 48109
//	            melsey@umich.edu
//
//              James A. Rossmanith
//              Department of Mathematics
//              University of Wisconsin
//              Madison, WI 53706
//              rossmani@math.wisc.edu
// =========================================================================

int main(int argc, char* argv[])
{
    //
    // NOTE: You should not have to modify this part of the code.
    //       To change parameters, modify the following files:
    //            1. input2D.data
    //            2. SignedDistance.cpp
    //            3. GridSpacing.cpp
    //
  
    // Get current time
    double time1 = time(NULL);

    // Call the mesh generator
    int meshgen2D(int argc,char**argv);
    int m = meshgen2D(argc, argv);
    
    // Get current time
    double time2 = time(NULL);

    // Output elapsed time
    cout << setprecision(5);
    cout << " Total elapsed time in seconds = " << setw(11) 
         << scientific << time2-time1 << endl << endl;

    return m;
}
