#include <string>
#include "meshdefs.h"

//
// Post-processing BEFORE all boundary and edge information is computed
//
void MeshPostProcess1(const string& GridType,double& h0, int& numpts, 
		      int& numtri, point*& p, triangle*& t, double*& area, 
		      double*& cdual, double (*SignedDistance)(point))
{

}
