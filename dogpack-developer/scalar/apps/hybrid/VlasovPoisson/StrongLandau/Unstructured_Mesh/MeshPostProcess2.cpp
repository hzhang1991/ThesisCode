#include "meshdefs.h"
#include "mesh.h"

//
// Post-processing AFTER all boundary and edge information is computed
// Here we modify "ghost_link" in order to enforce periodic boundary conditions
void MeshPostProcess2(char*& GridType, 
		      mesh& Mesh, 
		      double (*SignedDistance)(point))
{
  const int  NumPhysElems = Mesh.get_NumPhysElems();
  const int NumGhostElems = Mesh.get_NumGhostElems();
  void QuickSort(double*& a, int*& index, int lo, int hi);
  iTensor1 ghost_link(NumGhostElems);

  // ------------------------------------------------
  // Re-label nodes in each ghost element so that the
  // labeling is identical to the physical element
  // assiociated with current element.
  // ------------------------------------------------
  dTensor1 n1(2);
  dTensor1 n2(2);
  dTensor1 n3(2);
  dTensor1 m1(2);
  dTensor1 m2(2);
  dTensor1 m3(2);

  for (int i=1; i<=NumGhostElems; i++)
    {       
      int j = Mesh.get_ghost_link(i);
      int k = NumPhysElems+i;
      
      double xp1 = Mesh.get_node(Mesh.get_tnode(j,1),1);
      double yp1 = Mesh.get_node(Mesh.get_tnode(j,1),2);
      
      double xp2 = Mesh.get_node(Mesh.get_tnode(j,2),1);
      double yp2 = Mesh.get_node(Mesh.get_tnode(j,2),2);

      double xp3 = Mesh.get_node(Mesh.get_tnode(j,3),1);
      double yp3 = Mesh.get_node(Mesh.get_tnode(j,3),2);

      double L1 = sqrt(pow(yp3-yp2,2)+pow(xp2-xp3,2));
      n1.set(1, (yp3-yp2)/L1 );
      n1.set(2, (xp2-xp3)/L1 );
      
      double L2 = sqrt(pow(yp1-yp3,2)+pow(xp3-xp1,2));
      n2.set(1, (yp1-yp3)/L2 );
      n2.set(2, (xp3-xp1)/L2 );
      
      double L3 = sqrt(pow(yp2-yp1,2)+pow(xp1-xp2,2));
      n3.set(1, (yp2-yp1)/L3 );
      n3.set(2, (xp1-xp2)/L3 );

      double xg1 = Mesh.get_node(Mesh.get_tnode(k,1),1);
      double yg1 = Mesh.get_node(Mesh.get_tnode(k,1),2);
	    
      double xg2 = Mesh.get_node(Mesh.get_tnode(k,2),1);
      double yg2 = Mesh.get_node(Mesh.get_tnode(k,2),2);

      double xg3 = Mesh.get_node(Mesh.get_tnode(k,3),1);
      double yg3 = Mesh.get_node(Mesh.get_tnode(k,3),2);

      L1 = sqrt(pow(yg3-yg2,2)+pow(xg2-xg3,2));
      m1.set(1, (yg3-yg2)/L1 );
      m1.set(2, (xg2-xg3)/L1 );
      
      L2 = sqrt(pow(yg1-yg3,2)+pow(xg3-xg1,2));
      m2.set(1, (yg1-yg3)/L2 );
      m2.set(2, (xg3-xg1)/L2 );
      
      L3 = sqrt(pow(yg2-yg1,2)+pow(xg1-xg2,2));
      m3.set(1, (yg2-yg1)/L3 );
      m3.set(2, (xg1-xg2)/L3 );

      double dotprod1 = m1.get(1)*n1.get(1) + m1.get(2)*n1.get(2);
      double dotprod2 = m2.get(1)*n2.get(1) + m2.get(2)*n2.get(2);
      double dotprod3 = m3.get(1)*n3.get(1) + m3.get(2)*n3.get(2);

      if (fabs(dotprod1-1.0)>1.0e-12 || fabs(dotprod2-1.0)>1.0e-12 || fabs(dotprod3-1.0)>1.0e-12)
	{
	  xg1 = Mesh.get_node(Mesh.get_tnode(k,3),1);
	  yg1 = Mesh.get_node(Mesh.get_tnode(k,3),2);
	    
	  xg2 = Mesh.get_node(Mesh.get_tnode(k,1),1);
	  yg2 = Mesh.get_node(Mesh.get_tnode(k,1),2);

	  xg3 = Mesh.get_node(Mesh.get_tnode(k,2),1);
	  yg3 = Mesh.get_node(Mesh.get_tnode(k,2),2);

	  L1 = sqrt(pow(yg3-yg2,2)+pow(xg2-xg3,2));
	  m1.set(1, (yg3-yg2)/L1 );
	  m1.set(2, (xg2-xg3)/L1 );
	  
	  L2 = sqrt(pow(yg1-yg3,2)+pow(xg3-xg1,2));
	  m2.set(1, (yg1-yg3)/L2 );
	  m2.set(2, (xg3-xg1)/L2 );
	  
	  L3 = sqrt(pow(yg2-yg1,2)+pow(xg1-xg2,2));
	  m3.set(1, (yg2-yg1)/L3 );
	  m3.set(2, (xg1-xg2)/L3 );
	  
	  dotprod1 = m1.get(1)*n1.get(1) + m1.get(2)*n1.get(2);
	  dotprod2 = m2.get(1)*n2.get(1) + m2.get(2)*n2.get(2);
	  dotprod3 = m3.get(1)*n3.get(1) + m3.get(2)*n3.get(2);

	  if (fabs(dotprod1-1.0)<1.0e-12 && fabs(dotprod2-1.0)<1.0e-12 && fabs(dotprod3-1.0)<1.0e-12)
	    {
	      int tn1old = Mesh.get_tnode(k,1);
	      int tn2old = Mesh.get_tnode(k,2);
	      int tn3old = Mesh.get_tnode(k,3);

	      Mesh.set_tnode(k,1, tn3old );
	      Mesh.set_tnode(k,2, tn1old );
	      Mesh.set_tnode(k,3, tn2old );
	    }
	  else
	    {	      
	      xg1 = Mesh.get_node(Mesh.get_tnode(k,2),1);
	      yg1 = Mesh.get_node(Mesh.get_tnode(k,2),2);
	    
	      xg2 = Mesh.get_node(Mesh.get_tnode(k,3),1);
	      yg2 = Mesh.get_node(Mesh.get_tnode(k,3),2);
	      
	      xg3 = Mesh.get_node(Mesh.get_tnode(k,1),1);
	      yg3 = Mesh.get_node(Mesh.get_tnode(k,1),2);
	      
	      L1 = sqrt(pow(yg3-yg2,2)+pow(xg2-xg3,2));
	      m1.set(1, (yg3-yg2)/L1 );
	      m1.set(2, (xg2-xg3)/L1 );
	      
	      L2 = sqrt(pow(yg1-yg3,2)+pow(xg3-xg1,2));
	      m2.set(1, (yg1-yg3)/L2 );
	      m2.set(2, (xg3-xg1)/L2 );
	  
	      L3 = sqrt(pow(yg2-yg1,2)+pow(xg1-xg2,2));
	      m3.set(1, (yg2-yg1)/L3 );
	      m3.set(2, (xg1-xg2)/L3 );
	      
	      dotprod1 = m1.get(1)*n1.get(1) + m1.get(2)*n1.get(2);
	      dotprod2 = m2.get(1)*n2.get(1) + m2.get(2)*n2.get(2);
	      dotprod3 = m3.get(1)*n3.get(1) + m3.get(2)*n3.get(2);

	      if (fabs(dotprod1-1.0)<1.0e-12 && fabs(dotprod2-1.0)<1.0e-12 && fabs(dotprod3-1.0)<1.0e-12)
		{
		  int tn1old = Mesh.get_tnode(k,1);
		  int tn2old = Mesh.get_tnode(k,2);
		  int tn3old = Mesh.get_tnode(k,3);
		  
		  Mesh.set_tnode(k,1, tn2old );
		  Mesh.set_tnode(k,2, tn3old );
		  Mesh.set_tnode(k,3, tn1old );
		}
	      else
		{
		  printf("\n");
		  printf("ERROR in MeshPostProcess2.cpp: Should never get to this line.\n");
		  printf("\n");
		  exit(1);
		}
	    }
	}

      // Final check
      xg1 = Mesh.get_node(Mesh.get_tnode(k,1),1);
      yg1 = Mesh.get_node(Mesh.get_tnode(k,1),2);
	    
      xg2 = Mesh.get_node(Mesh.get_tnode(k,2),1);
      yg2 = Mesh.get_node(Mesh.get_tnode(k,2),2);

      xg3 = Mesh.get_node(Mesh.get_tnode(k,3),1);
      yg3 = Mesh.get_node(Mesh.get_tnode(k,3),2);

      L1 = sqrt(pow(yg3-yg2,2)+pow(xg2-xg3,2));
      m1.set(1, (yg3-yg2)/L1 );
      m1.set(2, (xg2-xg3)/L1 );
      
      L2 = sqrt(pow(yg1-yg3,2)+pow(xg3-xg1,2));
      m2.set(1, (yg1-yg3)/L2 );
      m2.set(2, (xg3-xg1)/L2 );
      
      L3 = sqrt(pow(yg2-yg1,2)+pow(xg1-xg2,2));
      m3.set(1, (yg2-yg1)/L3 );
      m3.set(2, (xg1-xg2)/L3 );

      dotprod1 = m1.get(1)*n1.get(1) + m1.get(2)*n1.get(2);
      dotprod2 = m2.get(1)*n2.get(1) + m2.get(2)*n2.get(2);
      dotprod3 = m3.get(1)*n3.get(1) + m3.get(2)*n3.get(2);

      if (fabs(dotprod1-1.0)>1.0e-12 || fabs(dotprod2-1.0)>1.0e-12 || fabs(dotprod3-1.0)>1.0e-12)
	{
	  printf("\n");
	  printf("ERROR in MeshPostProcess2.cpp: Should never get to this line.\n");
	  printf("  dotprod1 = %e\n",dotprod1);
	  printf("  dotprod2 = %e\n",dotprod2);
	  printf("  dotprod3 = %e\n",dotprod3);
	  printf("\n");
	  exit(1);
	}
    }
  // ------------------------------------------------
  
}
