#include "dogdefs.h"
#include "mesh.h"

// This is a user-supplied routine that sets the the boundary conditions
//
// This assumes perfect absorption (i.e. non-periodic!). 
// This routine is valid for the Vlasov beam problem only.
//
// For the Landau problems, this routine gets swapped out in the application
// Makefiles.
//
void SetBndValues_Unst(const mesh& Mesh, dTensor3* q, dTensor3* aux)
{   

    // problem information (If this were to be pulled from DogParams, these
    // numbers would NOT be correct!  The reason is that each quadrature point
    // was actually saved as a separate "equation")
    const int meqn = q->getsize(2);
    const int kmax = q->getsize(3);
    const int maux = aux->getsize(2);

    // Mesh information
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
        int j = Mesh.get_ghost_link(i);

        for (int m=1; m<=meqn; m++)
        for (int k=1; k<=kmax; k++)
        {
            q->set(i+NumPhysElems,m,k, 0.0 );
        }
    }

}
