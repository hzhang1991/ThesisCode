#include "mesh.h"
#include "FullMatrix.h"

void ConstructA_CG1(const mesh& Mesh, FullMatrix& A)
{
  const int NumPhysElems = Mesh.get_NumPhysElems();
  const int NumBndNodes  = Mesh.get_NumBndNodes();

  // Loop over all elements in the mesh
  for (int i=1; i<=NumPhysElems; i++)
    {
      // Information for element i
      iTensor1 tt(3);
      for (int k=1; k<=3; k++)
	{  tt.set(k, Mesh.get_tnode(i,k) );  }

      dTensor1 x(3);
      dTensor1 y(3);      
      for (int k=1; k<=3; k++)
	{  
	  x.set(k, Mesh.get_node(tt.get(k),1) );
	  y.set(k, Mesh.get_node(tt.get(k),2) );  
	}

      double Area = Mesh.get_area_prim(i);
      
      // Compute the three Lagrange polynomials on this element
      dTensor2 Mat(3,3);
      for (int k=1; k<=3; k++)
	{
	  Mat.set(k,1, 1.0  );
	  Mat.set(k,2, x.get(k) );
	  Mat.set(k,3, y.get(k) );
	}
      dTensor2 MatInv(3,3);
      void InvertM_CG1(const dTensor2& Mat,dTensor2& MatInv);
      InvertM_CG1(Mat,MatInv);

      // Compute the gradients of the three Lagrange polynomials on this element
      dTensor2 gphi(3,2);
      for (int k=1; k<=3; k++)
	{
	  gphi.set(k,1, MatInv.get(2,k) );
	  gphi.set(k,2, MatInv.get(3,k) );
	}

      // Entries of the stiffness matrix A
      for (int j=1; j<=3; j++)
	for (int k=1; k<=3; k++)
	  {
	    double tmp = A.get(tt.get(j),tt.get(k)) +
	      Area*(gphi.get(j,1)*gphi.get(k,1)+gphi.get(j,2)*gphi.get(k,2));
	    A.set(tt.get(j),tt.get(k), tmp );
	  }
    }

  // Replace boundary node equations by Dirichlet boundary condition enforcement
  for (int i=1; i<=NumBndNodes; i++)
    {
      const int j=Mesh.get_bnd_node(i);
      
      for (int k=1; k<=A.get_NumCols(); k++)
	{
	  A.set(j,k, 0.0 );	  
	}
      for (int k=1; k<=A.get_NumRows(); k++)
	{
	  A.set(k,j, 0.0 );
	}
      A.set(j,j, 1.0 );
    }

  // Get sparse structure representation
  A.Sparsify();

}


void InvertM_CG1(const dTensor2& Mat,
		 dTensor2& MatInv)
{
  const double x1 = Mat.get(1,2);
  const double y1 = Mat.get(1,3);
  const double x2 = Mat.get(2,2);
  const double y2 = Mat.get(2,3);
  const double x3 = Mat.get(3,2);
  const double y3 = Mat.get(3,3);
  const double detM = x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2;
  
  MatInv.set(1,1, (x2*y3-y2*x3)/detM );
  MatInv.set(1,2, (y1*x3-x1*y3)/detM );
  MatInv.set(1,3, (x1*y2-y1*x2)/detM );

  MatInv.set(2,1, (y2-y3)/detM );
  MatInv.set(2,2, (y3-y1)/detM );
  MatInv.set(2,3, (y1-y2)/detM );
  
  MatInv.set(3,1, (x3-x2)/detM );
  MatInv.set(3,2, (x1-x3)/detM );
  MatInv.set(3,3, (x2-x1)/detM );
}
