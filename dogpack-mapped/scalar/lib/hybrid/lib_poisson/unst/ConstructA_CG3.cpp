#include "mesh.h"
#include "FullMatrix.h"
#include "constants.h"

void ConstructA_CG3(const mesh& Mesh, FullMatrix& A)
{
  const int NumPhysElems = Mesh.get_NumPhysElems();
  const int NumBndNodes  = Mesh.get_SubNumBndNodes();
  const int Asize = Mesh.get_SubNumPhysNodes();

  assert_eq(Asize,A.get_NumRows());
  assert_eq(Asize,A.get_NumCols());
  
  dTensor1 A1(10);
  dTensor1 A2(10);
  dTensor1 A3(10);
  dTensor1 A4(10);
  dTensor1 A5(10);
  dTensor1 A6(10);
  dTensor1 A7(10);
  dTensor1 A8(10);
  dTensor1 A9(10);
  dTensor1 A10(10);

  A1.set(1,   0.0          );
  A1.set(2,   0.0          );
  A1.set(3,   0.0          );
  A1.set(4,   0.0          );
  A1.set(5,   0.0          );
  A1.set(6,   1.0          );
  A1.set(7,   0.0          );
  A1.set(8,   0.0          );
  A1.set(9,   0.0          );
  A1.set(10,  0.0          );
  
  A2.set(1,   0.5          );
  A2.set(2,  -1.5          );
  A2.set(3,   1.5          );
  A2.set(4,  -0.5          );
  A2.set(5,  -1.5          );
  A2.set(6,   0.0          );
  A2.set(7,   1.5          );
  A2.set(8,   0.0          );
  A2.set(9,   0.0          );
  A2.set(10,  0.0          );
      
  A3.set(1,   0.5          );
  A3.set(2,  -1.5          );
  A3.set(3,   0.0          );
  A3.set(4,   0.0          );
  A3.set(5,  -1.5          );
  A3.set(6,   0.0          );
  A3.set(7,   0.0          );
  A3.set(8,   1.5          );
  A3.set(9,   1.5          );
  A3.set(10, -0.5          );

  A4.set(1,   0.0          );
  A4.set(2,   4.5          );
  A4.set(3,  -4.5          );
  A4.set(4,   0.0          );
  A4.set(5,   4.5          );
  A4.set(6,  -9.0          );
  A4.set(7,   4.5          );
  A4.set(8,  -4.5          );
  A4.set(9,   4.5          );
  A4.set(10,  0.0          );
  
  A5.set(1,   0.0          );
  A5.set(2,   0.0          );
  A5.set(3,   0.0          );
  A5.set(4,   0.0          );
  A5.set(5,   4.5          );
  A5.set(6,  -9.0          );
  A5.set(7,   4.5          );
  A5.set(8,   0.0          );
  A5.set(9,   0.0          );
  A5.set(10,  0.0          );
  
  A6.set(1,   0.0          );
  A6.set(2,   4.5          );
  A6.set(3,   0.0          );
  A6.set(4,   0.0          );
  A6.set(5,   0.0          );
  A6.set(6,  -9.0          );
  A6.set(7,   0.0          );
  A6.set(8,   0.0          );
  A6.set(9,   4.5          );
  A6.set(10,  0.0          );
  
  A7.set(1,  -4.5          );
  A7.set(2,  13.5          );
  A7.set(3, -13.5          );
  A7.set(4,   4.5          );
  A7.set(5,   0.0          );
  A7.set(6,   0.0          );
  A7.set(7,   0.0          );
  A7.set(8,   0.0          );
  A7.set(9,   0.0          );
  A7.set(10,  0.0          );
  
  A8.set(1, -13.5          );
  A8.set(2,  27.0          );
  A8.set(3, -13.5          );
  A8.set(4,   0.0          );
  A8.set(5,  13.5          );
  A8.set(6, -27.0          );
  A8.set(7,  13.5          );
  A8.set(8,   0.0          );
  A8.set(9,   0.0          );
  A8.set(10,  0.0          );
  
  A9.set(1, -13.5          );
  A9.set(2,  13.5          );
  A9.set(3,   0.0          );
  A9.set(4,   0.0          );
  A9.set(5,  27.0          );
  A9.set(6, -27.0          );
  A9.set(7,   0.0          );
  A9.set(8, -13.5          );
  A9.set(9,  13.5          );
  A9.set(10,  0.0          );
  
  A10.set(1,  -4.5          );
  A10.set(2,   0.0          );
  A10.set(3,   0.0          );
  A10.set(4,   0.0          );
  A10.set(5,  13.5          );
  A10.set(6,   0.0          );
  A10.set(7,   0.0          );
  A10.set(8, -13.5          );
  A10.set(9,   0.0          );
  A10.set(10,  4.5          );

  dTensor2 spts(6,2);
  spts.set(1,1,  0.112615157582632 );
  spts.set(1,2,  0.112615157582632 );
  
  spts.set(2,1, -0.225230315165263 );
  spts.set(2,2,  0.112615157582632 );
  
  spts.set(3,1,  0.112615157582632 );
  spts.set(3,2, -0.225230315165263 );
  
  spts.set(4,1, -0.241757119823562 );
  spts.set(4,2, -0.241757119823562 );
  
  spts.set(5,1,  0.483514239647126 );
  spts.set(5,2, -0.241757119823562 );
  
  spts.set(6,1, -0.241757119823562 );
  spts.set(6,2,  0.483514239647126 );
    
  dTensor1 wgts(6);
  wgts.set(1, 0.1116907948390055 );
  wgts.set(2, 0.1116907948390055 );
  wgts.set(3, 0.1116907948390055 );
  wgts.set(4, 0.0549758718276610 );
  wgts.set(5, 0.0549758718276610 );
  wgts.set(6, 0.0549758718276610 );
  
  // Loop over all elements in the mesh
  for (int i=1; i<=NumPhysElems; i++)
    {
      // Information for element i
      iTensor1 tt(10);
      for (int k=1; k<=10; k++)
	{  tt.set(k, Mesh.get_node_subs(i,k) );  }
      
      // Evaluate gradients of the Lagrange polynomials on Gauss quadrature points      
      dTensor2 gpx(10,6);
      dTensor2 gpy(10,6);
      
      for (int m=1; m<=6; m++)
	{
	  double  xi = spts.get(m,1);
	  double eta = spts.get(m,2);
	  
	  for (int k=1; k<=10; k++)
	    {
	      double gp_xi  = A2.get(k) + A4.get(k)*eta + 2.0*A5.get(k)*xi
		+ 3.0*A7.get(k)*xi*xi + 2.0*A8.get(k)*xi*eta + A9.get(k)*eta*eta;
	      double gp_eta = A3.get(k) + A4.get(k)*xi + 2.0*A6.get(k)*eta 
		+ A8.get(k)*xi*xi + 2.0*A9.get(k)*xi*eta + 3.0*A10.get(k)*eta*eta;

	      gpx.set(k,m, Mesh.get_jmat(i,1,1)*gp_xi
		         + Mesh.get_jmat(i,1,2)*gp_eta );
	      gpy.set(k,m, Mesh.get_jmat(i,2,1)*gp_xi
		         + Mesh.get_jmat(i,2,2)*gp_eta );
	    }
	}

      // Entries of the stiffness matrix A
      double Area = Mesh.get_area_prim(i);
      for (int j=1; j<=10; j++)
	for (int k=1; k<=10; k++)
	  {
	    double tmp = A.get(tt.get(j),tt.get(k));
	    for (int m=1; m<=6; m++)
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
