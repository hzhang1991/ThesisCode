#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"

double GetCFL_Unst(double dt, const mesh& Mesh,
        const dTensor3& aux, const dTensor1& smax)
{

    double cfl=-100.0;
    int NumPhysElems = Mesh.get_NumPhysElems();
    int NumEdges = Mesh.get_NumEdges();

    for (int i=1; i<=NumPhysElems; i++)
    {
        double Area  = Mesh.get_area_prim(i);
        int edge1 = Mesh.get_tedge(i,1);
        int edge2 = Mesh.get_tedge(i,2);
        int edge3 = Mesh.get_tedge(i,3);

        double tmp = Max(Max(smax.get(edge1),smax.get(edge2)),
                smax.get(edge3));

        cfl = Max(0.5*dt*tmp/Area, cfl);
    }

    return cfl;

}
