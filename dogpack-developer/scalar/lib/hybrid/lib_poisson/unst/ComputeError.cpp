#include "mesh.h"
#include "tensors.h"

void ComputeError(const int space_order,
		  const mesh& Mesh,
		  const dTensor1& phi,
		  const dTensor2& E1,
		  const dTensor2& E2,
		  void (*PhiFunc)(const dTensor2& xpts,dTensor2& phi_ex),
		  void (*EfieldFunc)(const dTensor2& xpts,dTensor2& Efield_ex))
{
  
  // Potential
  const int NumPhysNodes = phi.getsize();//Mesh.get_NumPhysNodes();
  dTensor2 xpts(NumPhysNodes,2);
  dTensor2 phi_ex(NumPhysNodes,1);
  double phi_err;
  double phi_rel;
  
  switch(space_order)
    {
    case 2:
      for (int i=1; i<=NumPhysNodes; i++)
	{  
	  xpts.set(i,1, Mesh.get_node(i,1) );
	  xpts.set(i,2, Mesh.get_node(i,2) );
	}
      PhiFunc(xpts,phi_ex);
      
      phi_err = 0.0;
      phi_rel = 0.0;
      for (int i=1; i<=NumPhysNodes; i++)
	{ 
	  phi_rel = phi_rel + pow(phi_ex.get(i,1),2);
	  phi_err = phi_err + pow(phi_ex.get(i,1)-phi.get(i),2);
	}
      phi_err = sqrt(phi_err/phi_rel);
      break;

    case 3:
      for (int i=1; i<=NumPhysNodes; i++)
	{  
	  xpts.set(i,1, Mesh.get_sub_node(i,1) );
	  xpts.set(i,2, Mesh.get_sub_node(i,2) );
	}
      PhiFunc(xpts,phi_ex);
      
      phi_err = 0.0;
      phi_rel = 0.0;
      for (int i=1; i<=NumPhysNodes; i++)
	{ 
	  phi_rel = phi_rel + pow(phi_ex.get(i,1),2);
	  phi_err = phi_err + pow(phi_ex.get(i,1)-phi.get(i),2);
	}
      phi_err = sqrt(phi_err/phi_rel);
      break;
    }
  
  // Electric field components
  void L2Project_Unst(const int istart, 
		      const int iend, 
		      const int QuadOrder,		    
		      const int BasisOrder_fout,
		      const mesh& Mesh, 
		      dTensor3* fout, 
		      void (*Func)(const dTensor2&,dTensor2&));
  
  const int NumElems     = Mesh.get_NumElems();
  const int NumPhysElems = Mesh.get_NumPhysElems();
  const int kmax = E1.getsize(2);

  dTensor3 Efield_ex(NumElems,2,kmax);
  L2Project_Unst(1,NumElems,space_order,space_order,
		 Mesh,&Efield_ex,EfieldFunc);
  
  double E1_err = 0.0;
  double E1_rel = 0.0;
  double E2_err = 0.0;
  double E2_rel = 0.0;

  for (int i=1; i<=NumPhysElems; i++)
    {
      double Area = Mesh.get_area_prim(i);
      double tmp1 = 0.0;
      double tmp2 = 0.0;
      double tmp1_rel = 0.0;
      double tmp2_rel = 0.0;
      for (int k=1; k<=kmax; k++)
	{
	  tmp1 = tmp1 + pow((E1.get(i,k)-Efield_ex.get(i,1,k)),2);
	  tmp2 = tmp2 + pow((E2.get(i,k)-Efield_ex.get(i,2,k)),2);
	  tmp1_rel = tmp1_rel + pow((Efield_ex.get(i,1,k)),2);
	  tmp2_rel = tmp2_rel + pow((Efield_ex.get(i,2,k)),2);
	}
      E1_err = E1_err + Area*tmp1;
      E2_err = E2_err + Area*tmp2;
      E1_rel = E1_rel + Area*tmp1_rel;
      E2_rel = E2_rel + Area*tmp2_rel;
    }
  E1_err = sqrt(E1_err/E1_rel);
  E2_err = sqrt(E2_err/E2_rel);

  // Summary
  printf("  |----------------------------\n");
  printf("  | Errors:\n");
  printf("  |----------------------------\n");
  printf("  |  phi_err = %e\n",phi_err);
  printf("  |   E1_err = %e\n",E1_err);
  printf("  |   E2_err = %e\n",E2_err);
  printf("  |----------------------------\n");
  printf("\n");
}
