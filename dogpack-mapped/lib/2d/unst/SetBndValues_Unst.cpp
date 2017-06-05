#include "dogdefs.h"
#include "mesh.h"

// This is a user-supplied routine that sets the the boundary conditions
//
void SetBndValues_Unst(const mesh& Mesh, dTensor3* q, dTensor3* aux)
{   
  int meqn = q->getsize(2);
  int kmax = q->getsize(3);
  int maux = aux->getsize(2);
  int      NumElems = Mesh.get_NumElems();
  int  NumPhysElems = Mesh.get_NumPhysElems();
  int NumGhostElems = Mesh.get_NumGhostElems();
  int      NumNodes = Mesh.get_NumNodes();
  int  NumPhysNodes = Mesh.get_NumPhysNodes();
  int      NumEdges = Mesh.get_NumEdges();
        
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
	    q->set(i+NumPhysElems,m,k, q->get(j,m,k) );
	  }
      
      for (int m=1; m<=maux; m++)
	for (int k=1; k<=kmax; k++)
	  {
	    aux->set(i+NumPhysElems,m,k, aux->get(j,m,k) );
	  }
    }

}
