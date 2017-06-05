#include "mesh.h"
#include "constants.h"

void SetInterpNodes(const int space_order,
		    const mesh& Mesh,
		    const double hmin,
		    iTensor2& interp_nodes,
		    dTensor3& interp_mat)
{
  const int NumPhysElems = Mesh.get_NumPhysElems();
  const int NumPhysNodes = Mesh.get_NumPhysNodes();
  const int NumNodes     = Mesh.get_NumNodes();
  void QuickSort_int(iTensor1& a, iTensor1& index, int lo, int hi);
  void QuickSort_double(dTensor1& a, iTensor1& index, int lo, int hi);
  void GaussElimMatrixInv(dTensor2 inMat,dTensor2& outMat);
  
  // Find neighboring elements
  iTensor2 new_elem(NumPhysElems,4);
  void GetNeighboringElems(const double hmin, const mesh& Mesh, iTensor2& new_elem);
  GetNeighboringElems(hmin,Mesh,new_elem);

  switch(space_order)
    {
    case 2:
	  
      for (int i=1; i<=NumPhysElems; i++)
	{
	  // Grab all nodes ...
	  iTensor1 n(12);
	  iTensor1 nindex(12);
	  int num = 0;
	  for (int j=1; j<=4; j++)
	    {
	      for (int k=1; k<=3; k++)
		{
		  num = num+1;
		  n.set(num, Mesh.get_tnode(new_elem.get(i,j),k) );
		  nindex.set(num, num );
		}
	    }

	  // ... then sort all grabbed nodes ...
	  QuickSort_int(n,nindex,1,12);

	  // ... then only pick unique ones and store int "interp_nodes"
	  interp_nodes.set(i,1, n.get(1) );
	  int ijk=1;
	  for (int k=2; k<=12; k++)
	    {
	      int mtmp1 = interp_nodes.get(i,ijk);
	      int mtmp2 = n.get(k);
	      if (mtmp1!=mtmp2)
		{
		  ijk = ijk+1;
		  interp_nodes.set(i,ijk, mtmp2 );
		}
	    }
	  if (ijk!=6)
	    {
	      printf("\n");
	      printf(" Error 3 \n");
	      printf("\n");
	      exit(1);
	    }

	  /*
	  for (int j=1; j<=6; j++)
	    {
	      printf(" interp_nodes.get(%i,%i) = %i\n",i,j,interp_nodes.get(i,j) );
	    }
	  printf("\n");
	  */

	  // Compute inversion of interpolation matrix and store
	  // in the matrix "interp_mat"	  
	  dTensor2 Mat(6,6);
	  dTensor2 MatInv(6,6);
	  
          for (int k=1; k<=6; k++)
            {
	      int nnn = interp_nodes.get(i,k);
	      double x = Mesh.get_node(nnn,1);
	      double y = Mesh.get_node(nnn,2);

              Mat.set(k,1, 1.0 );
              Mat.set(k,2, x   );
              Mat.set(k,3, y   );
              Mat.set(k,4, x*y );
              Mat.set(k,5, x*x );
              Mat.set(k,6, y*y );
            }	  
	    
	  GaussElimMatrixInv(Mat,MatInv);
	  for (int j=1; j<=6; j++)
	    for (int k=1; k<=6; k++)
	      {
		interp_mat.set(i,j,k, MatInv.get(j,k) );
	      }	  
	  		   
	}
      break;

    case 3:

      for (int i=1; i<=NumPhysElems; i++)
	{
	  // Grab all nodes ...
	  iTensor1 n(24);
	  iTensor1 nindex(24);
	  int num = 0;
	  for (int j=1; j<=4; j++)
	    {
	      for (int k=1; k<=6; k++)
		{
		  num = num+1;
		  n.set(num, Mesh.get_node_subs(new_elem.get(i,j),k) );
		  nindex.set(num, num );
		}
	    }

	  // ... then sort all grabbed nodes ...
	  QuickSort_int(n,nindex,1,24);

	  // ... then only pick unique ones and store int "interp_nodes"
	  interp_nodes.set(i,1, n.get(1) );
	  int ijk=1;
	  for (int k=2; k<=24; k++)
	    {
	      int mtmp1 = interp_nodes.get(i,ijk);
	      int mtmp2 = n.get(k);
	      if (mtmp1!=mtmp2)
		{
		  ijk = ijk+1;
		  interp_nodes.set(i,ijk, mtmp2 );
		}
	    }
	  if (ijk!=15)
	    {
	      printf("\n");
	      printf(" Error 3 \n");
	      printf("\n");
	      exit(1);
	    }

	  /*
	  for (int j=1; j<=6; j++)
	    {
	      printf(" interp_nodes.get(%i,%i) = %i\n",i,j,interp_nodes.get(i,j) );
	    }
	  printf("\n");
	  */

	  // Compute inversion of interpolation matrix and store
	  // in the matrix "interp_mat"	  
	  dTensor2 Mat(15,15);
	  dTensor2 MatInv(15,15);
	  
          for (int k=1; k<=15; k++)
            {
	      int nnn = interp_nodes.get(i,k);
	      double x = Mesh.get_sub_node(nnn,1);
	      double y = Mesh.get_sub_node(nnn,2);
	      double x2 = x*x;
	      double x3 = x*x2;
	      double x4 = x*x3;
	      double y2 = y*y;
	      double y3 = y*y2;
	      double y4 = y*y3;
	      
              Mat.set(k,1,  1.0   );
              Mat.set(k,2,  x     );
              Mat.set(k,3,  y     );
              Mat.set(k,4,  x*y   );
              Mat.set(k,5,  x2    );
              Mat.set(k,6,  y2    );
	      Mat.set(k,7,  x3    );
	      Mat.set(k,8,  x2*y  );
	      Mat.set(k,9,  x*y2  );
	      Mat.set(k,10, y3    );
	      Mat.set(k,11, x4    );
	      Mat.set(k,12, x3*y  );
	      Mat.set(k,13, x2*y2 );
	      Mat.set(k,14, x*y3  );
	      Mat.set(k,15, y4    );
            }	  

	  /*
	  printf("\n");
	  printf("  M := diag(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0): \n");
	  for (int j=1; j<=15; j++)
	    for (int k=1; k<=15; k++)
	      {
		printf(" M[%i,%i]:= %e;\n",j,k,Mat.get(j,k));
	      }
	  printf("\n");
	  */
	  GaussElimMatrixInv(Mat,MatInv);
	  for (int j=1; j<=15; j++)
	    for (int k=1; k<=15; k++)
	      {
		interp_mat.set(i,j,k, MatInv.get(j,k) );
	      }	  
	  		   
	}      
      break;
      
    default:
      printf("\n");
      printf(" Error in SetInterpNodes.cpp: space_order value is not implemented ...\n");
      printf("     space_order = %i\n",space_order);
      printf("\n");
      exit(1);
      break;
    }

}


void GetNeighboringElems(const double hmin, const mesh& Mesh, iTensor2& new_elem)
{
  const int NumPhysElems = Mesh.get_NumPhysElems();
  void QuickSort_int(iTensor1& a, iTensor1& index, int lo, int hi);
  void QuickSort_double(dTensor1& a, iTensor1& index, int lo, int hi);

  for (int i=1; i<=NumPhysElems; i++)
    {	  
      // Find 3 neighboring elements
      new_elem.set(i,1, 0 );
      new_elem.set(i,2, 0 );
      new_elem.set(i,3, 0 );
      new_elem.set(i,4, i );
      int num_outside = 0;
      int num_inside = 0;
      iTensor1 outside(3);
      iTensor1 inside(3);
      dTensor2 xyc(4,2);
      
      int mn1=Mesh.get_tnode(i,1);
      int mn2=Mesh.get_tnode(i,2);
      int mn3=Mesh.get_tnode(i,3);
      xyc.set(1,1, onethird*(Mesh.get_node(mn1,1) +
			     Mesh.get_node(mn2,1) +
			     Mesh.get_node(mn3,1) ) );
      xyc.set(1,2, onethird*(Mesh.get_node(mn1,2) +
			     Mesh.get_node(mn2,2) +
			     Mesh.get_node(mn3,2) ) );
      
      for (int k=1; k<=3; k++)
	{
	  int e = Mesh.get_tedge(i,k);
	  
	  int tmp1 = Mesh.get_eelem(e,1);
	  int tmp2 = Mesh.get_eelem(e,2);
	  
	  if (tmp1==i)
	    {  new_elem.set(i,k, tmp2 );  }
	  else
	    {  new_elem.set(i,k, tmp1 );  }
	  
	  mn1=Mesh.get_tnode(new_elem.get(i,k),1);
	  mn2=Mesh.get_tnode(new_elem.get(i,k),2);
	  mn3=Mesh.get_tnode(new_elem.get(i,k),3);
	  xyc.set(k+1,1, onethird*(Mesh.get_node(mn1,1) +
				   Mesh.get_node(mn2,1) +
				   Mesh.get_node(mn3,1) ) );
	  xyc.set(k+1,2, onethird*(Mesh.get_node(mn1,2) +
				   Mesh.get_node(mn2,2) +
				   Mesh.get_node(mn3,2) ) ); 
	  
	  if (new_elem.get(i,k) > NumPhysElems)
	    {
	      num_outside = num_outside + 1;
	      outside.set(num_outside, k );
	    }
	  else
	    {
	      num_inside = num_inside + 1;
	      inside.set(num_inside, k );
	    }
	}
      
      if ((num_inside+num_outside)!=3)
	{
	  printf("\n");
	  printf(" ERROR in SetInterpNodes: \n");
	  printf("    num_outside = %i\n",num_outside);
	  printf("     num_inside = %i\n",num_inside);
	  printf("              i = %i\n",i);
	  printf("\n");
	  exit(1);
	}
      
      /*
      // Scenarios:
      //   1. All elements inside computational domain
      //   2. One element is outside computational domain
      //   3. Two elements are outside computational domain
      if( num_outside==1 )
	{
	  int tprob1 = new_elem.get(i,outside.get(1));
	  int tfine1 = new_elem.get(i,inside.get(1));
	  int tfine2 = new_elem.get(i,inside.get(2));
	  
	  int mfound = 0;
	  int mmm = 0;
	  int mtmp;
	  while(mfound==0)
	    {
	      mmm=mmm+1;
	      if (mmm>3)		    
		{
		  printf("\n");
		  printf(" Error 1 \n");
		  printf("\n");
		  exit(1);
		}
	      
	      int e = Mesh.get_tedge(tfine1,mmm);
	      
	      int tmp1 = Mesh.get_eelem(e,1);
	      int tmp2 = Mesh.get_eelem(e,2);
	      
	      if (tmp1==tfine1)
		{  mtmp = tmp2;  }
	      else
		{  mtmp = tmp1;  }
	      
	      if ( (mtmp<=NumPhysElems) && (mtmp!=i) )
		{
		  int mn1=Mesh.get_tnode(mtmp,1);
		  int mn2=Mesh.get_tnode(mtmp,2);
		  int mn3=Mesh.get_tnode(mtmp,3);
		  xyc.set(outside.get(1)+1,1, onethird*(Mesh.get_node(mn1,1) +
							Mesh.get_node(mn2,1) +
							Mesh.get_node(mn3,1) ) );
		  xyc.set(outside.get(1)+1,2, onethird*(Mesh.get_node(mn1,2) +
							Mesh.get_node(mn2,2) +
							Mesh.get_node(mn3,2) ) ); 
		  
		  dTensor1 xsort(4);
		  dTensor1 ysort(4);
		  iTensor1 xindex(4);
		  iTensor1 yindex(4);
		  for (int ijk=1; ijk<=4; ijk++)
		    {
		      xsort.set(ijk, xyc.get(ijk,1) );
		      xindex.set(ijk, ijk );
		      ysort.set(ijk, xyc.get(ijk,2) );
		      yindex.set(ijk, ijk );
		    }
		  
		  QuickSort_double(xsort,xindex,1,4);		      
		  QuickSort_double(ysort,yindex,1,4);
		  
		  int xunique = 1;
		  int yunique = 1;
		  for (int ijk=2; ijk<=4; ijk++)
		    {
		      double a1 = xsort.get(ijk-1);
		      double a2 = xsort.get(ijk);
		      
		      if (fabs(a1-a2)>hmin)
			{
			  xunique = xunique + 1;
			}
		      
		      a1 = ysort.get(ijk-1);
		      a2 = ysort.get(ijk);
		      
		      if (fabs(a1-a2)>hmin)
			{
			  yunique = yunique + 1;
			}
		    }
		  
		  //printf(" i = %i, xunique = %i,  yunique = %i \n",i,xunique,yunique);
		  //printf(" xsort = %e, %e, %e, %e\n",xsort.get(1),xsort.get(2),xsort.get(3),xsort.get(4));
		  //printf(" ysort = %e, %e, %e, %e\n",ysort.get(1),ysort.get(2),ysort.get(3),ysort.get(4));
		  
		  if (xunique>2 && yunique>2)
		    {
		      new_elem.set(i,outside.get(1), mtmp );
		      mfound = 1;
		    }
		}
	      
	      if (mfound==0)
		{
		  e = Mesh.get_tedge(tfine2,mmm);
		  
		  tmp1 = Mesh.get_eelem(e,1);
		  tmp2 = Mesh.get_eelem(e,2);
		      
		  if (tmp1==tfine2)
		    {  mtmp = tmp2;  }
		  else
		    {  mtmp = tmp1;  }
		  
		  if ( (mtmp<=NumPhysElems) && (mtmp!=i) )
		    {
		      int mn1=Mesh.get_tnode(mtmp,1);
		      int mn2=Mesh.get_tnode(mtmp,2);
		      int mn3=Mesh.get_tnode(mtmp,3);
		      xyc.set(outside.get(1)+1,1, onethird*(Mesh.get_node(mn1,1) +
							    Mesh.get_node(mn2,1) +
							    Mesh.get_node(mn3,1) ) );
		      xyc.set(outside.get(1)+1,2, onethird*(Mesh.get_node(mn1,2) +
							    Mesh.get_node(mn2,2) +
							    Mesh.get_node(mn3,2) ) ); 
		      
		      dTensor1 xsort(4);
		      dTensor1 ysort(4);
		      iTensor1 xindex(4);
		      iTensor1 yindex(4);
		      for (int ijk=1; ijk<=4; ijk++)
			{
			  xsort.set(ijk, xyc.get(ijk,1) );
			  xindex.set(ijk, ijk );
			  ysort.set(ijk, xyc.get(ijk,2) );
			  yindex.set(ijk, ijk );
			}
		      
		      QuickSort_double(xsort,xindex,1,4);		      
		      QuickSort_double(ysort,yindex,1,4);
		      
		      int xunique = 1;
		      int yunique = 1;
		      for (int ijk=2; ijk<=4; ijk++)
			{
			  double a1 = xsort.get(ijk-1);
			  double a2 = xsort.get(ijk);
			      
			  if (fabs(a1-a2)>hmin)
			    {
			      xunique = xunique + 1;
			    }
			  
			  a1 = ysort.get(ijk-1);
			  a2 = ysort.get(ijk);
			  
			  if (fabs(a1-a2)>hmin)
			    {
			      yunique = yunique + 1;
			    }
			}
		      
		      //printf(" i = %i, xunique = %i,  yunique = %i \n",i,xunique,yunique);
		      //printf(" xsort = %e, %e, %e, %e\n",xsort.get(1),xsort.get(2),xsort.get(3),xsort.get(4));
		      //printf(" ysort = %e, %e, %e, %e\n",ysort.get(1),ysort.get(2),ysort.get(3),ysort.get(4));  
		      
		      if (xunique>2 && yunique>2)
			{
			  new_elem.set(i,outside.get(1), mtmp );
			  mfound = 1;
			}
		    }
		}		  
	    }
	}
      else if (num_outside==2 )
	{
	  int tprob1 = new_elem.get(i,outside.get(1));
	  int tprob2 = new_elem.get(i,outside.get(2));
	  int tfine = new_elem.get(i,inside.get(1));
	  
	  int e = Mesh.get_tedge(i,inside.get(1));
	  
	  int tmp1 = Mesh.get_eelem(e,1);
	  int tmp2 = Mesh.get_eelem(e,2);
	  
	  int etmp1 = Mesh.get_tedge(tfine,1);
	  int etmp2 = Mesh.get_tedge(tfine,2);
	  int etmp3 = Mesh.get_tedge(tfine,3);
	  
	  int e1;
	  int e2;
	  
	  if (etmp1==e)
	    {
	      e1 = etmp2;
	      e2 = etmp3;
	    }
	  else if (etmp2==e)
	    {
	      e1 = etmp1;
	      e2 = etmp3;
	    }
	  else if (etmp3==e)
	    {
	      e1 = etmp1;
	      e2 = etmp2;
	    }
	  
	  tmp1 = Mesh.get_eelem(e1,1);
	  tmp2 = Mesh.get_eelem(e1,2);

	  if (tmp1==tfine)
	    {
	      new_elem.set(i,outside.get(1), tmp2 );
	    }
	  else 
	    {
	      new_elem.set(i,outside.get(1), tmp1 );
	    }
	  
	  tmp1 = Mesh.get_eelem(e2,1);
	  tmp2 = Mesh.get_eelem(e2,2);
	  
	  if (tmp1==tfine)
	    {
	      new_elem.set(i,outside.get(2), tmp2 );
	    }
	  else 
	    {
	      new_elem.set(i,outside.get(2), tmp1 );
	    }
	  
	}
      */
      
      iTensor1 new_elem_tmp(4);
      iTensor1 new_elem_index(4);
      for (int m=1; m<=4; m++)
	{  
	  new_elem_index.set(m, m );  
	  new_elem_tmp.set(m, new_elem.get(i,m) );
	}
      QuickSort_int(new_elem_tmp,new_elem_index,1,4);
      for (int m=1; m<=4; m++)
	{
	  new_elem.set(i,m, new_elem_tmp.get(m) );
	}

      /*
      if (num_outside>0)
	{
	  for (int m=1; m<=4; m++)
	    {
	      printf("(%i,%i) = %i    ",i,m,new_elem.get(i,m));
	    }
	  printf("\n");
	}
      */
      
      for (int m=1; m<=3; m++)
	{
	  int tmp1 = new_elem.get(i,m);
	  int tmp2 = new_elem.get(i,m+1);
	  
	  if (tmp1==tmp2)
	    { 		 
	      printf("\n");
	      printf(" Error 2 \n");
	      printf("\n");
	      exit(1);		
	    }
	}
      
      /*
	for (int j=1; j<=4; j++)
	{
	printf(" i = %i, new_elem.get(%i) = %i\n",i,j,new_elem.get(j) );
	}
	printf("\n");
      */
    }

}
