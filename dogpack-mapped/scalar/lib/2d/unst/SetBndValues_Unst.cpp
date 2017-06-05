#include "dogdefs.h"
#include "mesh.h"

// This is a user-supplied routine that sets the the boundary conditions
//
void SetBndValues_Unst(const mesh& Mesh, dTensor3* q, dTensor3* aux)
{   

    const int meqn = q->getsize(2);
    const int kmax = q->getsize(3);
    const int maux = aux->getsize(2);
    const int      NumElems = Mesh.get_NumElems();
    const int  NumPhysElems = Mesh.get_NumPhysElems();
    const int NumGhostElems = Mesh.get_NumGhostElems();
    const int      NumNodes = Mesh.get_NumNodes();
    const int  NumPhysNodes = Mesh.get_NumPhysNodes();
    const int      NumEdges = Mesh.get_NumEdges();

    // ----------------------------------------
    // Loop over each ghost cell element and
    // place the correct information into
    // these elements
    // ----------------------------------------
    for (int i=1; i<=NumGhostElems; i++)
    {
        const int j = Mesh.get_ghost_link(i);
        for (int m=1; m<=meqn; m++)
            for (int k=1; k<=kmax; k++)
            {
                q->set(i+NumPhysElems,m,k, q->get(j,m,k) );
            }

        for (int m=1; m<=maux; m++)
            for (int k=1; k<=kmax; k++)
            {
                aux->set(i+NumPhysElems,m,k, aux->get(j,m,k) );
            }
    }

}
