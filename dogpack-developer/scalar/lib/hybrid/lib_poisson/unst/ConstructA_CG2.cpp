#include "mesh.h"
#include "FullMatrix.h"
#include "constants.h"

void ConstructA_CG2(const mesh& Mesh, FullMatrix& A)
{
  const int NumPhysElems = Mesh.get_NumPhysElems();
  const int NumBndNodes  = Mesh.get_SubNumBndNodes();
  const int Asize = Mesh.get_SubNumPhysNodes();

  assert_eq(Asize,A.get_NumRows());
  assert_eq(Asize,A.get_NumCols());
  
  dTensor1 A1(6);
  dTensor1 A2(6);
  dTensor1 A3(6);
  dTensor1 A4(6);
  dTensor1 A5(6);
  dTensor1 A6(6);

  A1.set(1, -oneninth     );
  A1.set(2,  4.0*oneninth );
  A1.set(3, -oneninth     );
  A1.set(4,  4.0*oneninth );
  A1.set(5,  4.0*oneninth );
  A1.set(6, -oneninth     );
  
  A2.set(1, -onethird     );
  A2.set(2,  0.0          );
  A2.set(3,  onethird     );
  A2.set(4, -4.0*onethird );
  A2.set(5,  4.0*onethird );
  A2.set(6,  0.0          );
  
  A3.set(1, -onethird     );
  A3.set(2, -4.0*onethird );
  A3.set(3,  0.0          );
  A3.set(4,  0.0          );
  A3.set(5,  4.0*onethird );
  A3.set(6,  onethird     );
  
  A4.set(1,  4.0          );
  A4.set(2, -4.0          );
  A4.set(3,  0.0          );
  A4.set(4, -4.0          );
  A4.set(5,  4.0          );
  A4.set(6,  0.0          );

  A5.set(1,  2.0          );
  A5.set(2, -4.0          );
  A5.set(3,  2.0          );
  A5.set(4,  0.0          );
  A5.set(5,  0.0          );
  A5.set(6,  0.0          );
  
  A6.set(1,  2.0          );
  A6.set(2,  0.0          );
  A6.set(3,  0.0          );
  A6.set(4, -4.0          );
  A6.set(5,  0.0          );
  A6.set(6,  2.0          );

  dTensor2 spts(3,2);
  spts.set(1,1,  1.0/3.0 );
  spts.set(1,2, -1.0/6.0 );
  
  spts.set(2,1, -1.0/6.0 );
  spts.set(2,2, -1.0/6.0 );
  
  spts.set(3,1, -1.0/6.0 );
  spts.set(3,2,  1.0/3.0 );
  
  dTensor1 wgts(3);
  wgts.set(1, 1.0/6.0 );
  wgts.set(2, 1.0/6.0 );
  wgts.set(3, 1.0/6.0 );
  
  // Loop over all elements in the mesh
  for (int i=1; i<=NumPhysElems; i++)
    {
      // Information for element i
      iTensor1 tt(6);
      for (int k=1; k<=6; k++)
	{  tt.set(k, Mesh.get_node_subs(i,k) );  }
      
      // Evaluate gradients of the Lagrange polynomials on Gauss quadrature points      
      dTensor2 gpx(6,3);
      dTensor2 gpy(6,3);
      
      for (int m=1; m<=3; m++)
	{
	  double  xi = spts.get(m,1);
	  double eta = spts.get(m,2);
	  
	  for (int k=1; k<=6; k++)
	    {
	      double gp_xi  = A2.get(k) + 2.0*A5.get(k)*xi + A4.get(k)*eta;
	      double gp_eta = A3.get(k) + A4.get(k)*xi + 2.0*A6.get(k)*eta;

	      gpx.set(k,m, Mesh.get_jmat(i,1,1)*gp_xi
		         + Mesh.get_jmat(i,1,2)*gp_eta );
	      gpy.set(k,m, Mesh.get_jmat(i,2,1)*gp_xi
		         + Mesh.get_jmat(i,2,2)*gp_eta );
	    }
	}

      // Entries of the stiffness matrix A
      double Area = Mesh.get_area_prim(i);
      for (int j=1; j<=6; j++)
	for (int k=1; k<=6; k++)
	  {
	    double tmp = A.get(tt.get(j),tt.get(k));
	    for (int m=1; m<=3; m++)
	      {
		tmp = tmp + 2.0*Area*wgts.get(m)*(gpx.get(j,m)*gpx.get(k,m)+gpy.get(j,m)*gpy.get(k,m));
	      }
	    A.set(tt.get(j),tt.get(k), tmp );
	  }
    }

  // Replace boundary node equations by Dirichlet boundary condition enforcement
  for (int i=1; i<=NumBndNodes; i++)
    {
      const int j=Mesh.get_sub_bnd_node(i);
      
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
