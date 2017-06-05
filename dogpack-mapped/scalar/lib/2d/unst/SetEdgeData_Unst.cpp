#include "dogdefs.h"
#include "mesh.h"
#include "MonomialsToLegendre.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include <cmath>

void SetEdgeData_Unst(const mesh& Mesh, 
		      int NumQuadPoints, 
		      int NumBasisOrder, 
		      edge_data_Unst& EdgeData)
{
  // Quick error check
  if (NumQuadPoints<1 || NumQuadPoints>5 || NumBasisOrder<1 || NumBasisOrder>5)
    {
      printf(" \n");
      printf(" Error in SetEdgeData_Unst.cpp \n");
      printf("   NumQuadPoints and NumBasisOrder must be 1,2,3,4, or 5.\n");
      printf("     NumQuadPoints = %i\n",NumQuadPoints);
      printf("     NumBasisOrder = %i\n",NumBasisOrder);
      printf("\n");
      exit(1);
    }

  // ---------------------------------
  // Set quadrature weights and points
  // ---------------------------------
  switch( NumQuadPoints )
    {
    case 1:
      EdgeData.wgts1d->set(1, 2.0 );
      
      EdgeData.xpts1d->set(1, 0.0 );
      break;
      
    case 2:
      EdgeData.wgts1d->set(1,  1.0 );
      EdgeData.wgts1d->set(2,  1.0 );
      
      EdgeData.xpts1d->set(1,  1.0/sq3 );
      EdgeData.xpts1d->set(2, -1.0/sq3 );
      break;
      
    case 3:
      EdgeData.wgts1d->set(1,  5.0/9.0 );
      EdgeData.wgts1d->set(2,  8.0/9.0 );
      EdgeData.wgts1d->set(3,  5.0/9.0 );
      
      EdgeData.xpts1d->set(1,  sq3/sq5 );
      EdgeData.xpts1d->set(2,  0.0 );
      EdgeData.xpts1d->set(3, -sq3/sq5 );
      break;
      
    case 4:
      EdgeData.wgts1d->set(1, (18.0 - sq3*sq10)/36.0 );
      EdgeData.wgts1d->set(2, (18.0 + sq3*sq10)/36.0 );
      EdgeData.wgts1d->set(3, EdgeData.wgts1d->get(2) );
      EdgeData.wgts1d->set(4, EdgeData.wgts1d->get(1) );
      
      EdgeData.xpts1d->set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
      EdgeData.xpts1d->set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
      EdgeData.xpts1d->set(3, -EdgeData.xpts1d->get(2) );
      EdgeData.xpts1d->set(4, -EdgeData.xpts1d->get(1) );           
      break;
      
    case 5:      
      EdgeData.wgts1d->set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
      EdgeData.wgts1d->set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
      EdgeData.wgts1d->set(3, 128.0/225.0 );
      EdgeData.wgts1d->set(4, EdgeData.wgts1d->get(2) );
      EdgeData.wgts1d->set(5, EdgeData.wgts1d->get(1) );
      
      EdgeData.xpts1d->set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
      EdgeData.xpts1d->set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
      EdgeData.xpts1d->set(3,  0.0 );
      EdgeData.xpts1d->set(4, -EdgeData.xpts1d->get(2) );
      EdgeData.xpts1d->set(5, -EdgeData.xpts1d->get(1) );
      break;
    }

  // ---------------------------------
  // Legendre basis functions on the 
  // left and right of each edge
  // ---------------------------------
  const int NumEdges = Mesh.get_NumEdges();
  const int NumBasisComps = (NumBasisOrder*(NumBasisOrder+1))/2;
  dTensor1 xp1(3);
  dTensor1 yp1(3);
  dTensor1 xp2(3);
  dTensor1 yp2(3);
  dTensor1 xy1(2);
  dTensor1 xy2(2);
  dTensor1 mu1(NumBasisComps);
  dTensor1 mu2(NumBasisComps);

  for (int i=1; i<=NumEdges; i++)
    {   
      // Get edge information
      const double x1 = Mesh.get_edge(i,1);
      const double y1 = Mesh.get_edge(i,2);
      const double x2 = Mesh.get_edge(i,3);
      const double y2 = Mesh.get_edge(i,4);
      
      const int e1 = Mesh.get_eelem(i,1);
      const int e2 = Mesh.get_eelem(i,2);

      // Get element information about
      // the two elements that meet at
      // the current edge
      const double Area1 = Mesh.get_area_prim(e1);
      const double Area2 = Mesh.get_area_prim(e2);

      for (int k=1; k<=3; k++)
	{
	  xp1.set(k, Mesh.get_node(Mesh.get_tnode(e1,k),1) );
	  yp1.set(k, Mesh.get_node(Mesh.get_tnode(e1,k),2) );

	  xp2.set(k, Mesh.get_node(Mesh.get_tnode(e2,k),1) );
	  yp2.set(k, Mesh.get_node(Mesh.get_tnode(e2,k),2) );
	}

      const double xc1 = (xp1.get(1) + xp1.get(2) + xp1.get(3))/3.0;
      const double yc1 = (yp1.get(1) + yp1.get(2) + yp1.get(3))/3.0;
      const double xc2 = (xp2.get(1) + xp2.get(2) + xp2.get(3))/3.0;
      const double yc2 = (yp2.get(1) + yp2.get(2) + yp2.get(3))/3.0;

      // quadrature points on the edge
      for (int m=1; m<=NumQuadPoints; m++)
	{
	  // Take integration point s (in [-1,1])
	  // and map to physical domain
	  const double s = EdgeData.xpts1d->get(m);
	  const double x = x1 + 0.5*(s+1.0)*(x2-x1);
	  const double y = y1 + 0.5*(s+1.0)*(y2-y1);

	  // Take physical point (x,y)
	  // and map into the coordinates
	  // of the two triangles that are
	  // adjacent to the current edge
	  xy1.set(1, ((yp1.get(3)-yp1.get(1))*(x-xc1) 
		    + (xp1.get(1)-xp1.get(3))*(y-yc1))/(2.0*Area1) );
	  xy1.set(2, ((yp1.get(1)-yp1.get(2))*(x-xc1) 
		    + (xp1.get(2)-xp1.get(1))*(y-yc1))/(2.0*Area1) );
	  
	  xy2.set(1, ((yp2.get(3)-yp2.get(1))*(x-xc2) 
		    + (xp2.get(1)-xp2.get(3))*(y-yc2))/(2.0*Area2) );
	  xy2.set(2, ((yp2.get(1)-yp2.get(2))*(x-xc2) 
		    + (xp2.get(2)-xp2.get(1))*(y-yc2))/(2.0*Area2) );

	  // Evaluate monomials at locations xy1
	  double xi = xy1.get(1);
	  double xi2 = xi*xi;
	  double xi3 = xi*xi2;
	  double xi4 = xi*xi3;

	  double eta = xy1.get(2);
	  double eta2 = eta*eta;
	  double eta3 = eta*eta2;
	  double eta4 = eta*eta3;

	  switch( NumBasisOrder )
	    {
	    case 5:  // fifth order		    		    
	      mu1.set(15, eta4     );
	      mu1.set(14, xi4      );
	      mu1.set(13, xi2*eta2 );
	      mu1.set(12, eta3*xi  );
	      mu1.set(11, xi3*eta  );
	      
	    case 4:  // fourth order
	      mu1.set(10, eta3     );
	      mu1.set(9,  xi3      );
	      mu1.set(8,  xi*eta2  );
	      mu1.set(7,  eta*xi2  );
	      
	    case 3:  // third order
	      mu1.set(6,  eta2     );
	      mu1.set(5,  xi2      );
	      mu1.set(4,  xi*eta   );		    
	      
	    case 2:  // second order		    
	      mu1.set(3, eta       );
	      mu1.set(2, xi        );
	      
	    case 1:  // first order
	      mu1.set(1, 1.0       );
	      
	      break;		    
	    }
	  
	  // Evaluate monomials at locations xy2
	  xi = xy2.get(1);
	  xi2 = xi*xi;
	  xi3 = xi*xi2;
	  xi4 = xi*xi3;
	  
	  eta = xy2.get(2);
	  eta2 = eta*eta;
	  eta3 = eta*eta2;
	  eta4 = eta*eta3;

	  switch( NumBasisOrder )
	    {
	    case 5:  // fifth order		    		    
	      mu2.set(15, eta4     );
	      mu2.set(14, xi4      );
	      mu2.set(13, xi2*eta2 );
	      mu2.set(12, eta3*xi  );
	      mu2.set(11, xi3*eta  );
	      
	    case 4:  // fourth order
	      mu2.set(10, eta3     );
	      mu2.set(9,  xi3      );
	      mu2.set(8,  xi*eta2  );
	      mu2.set(7,  eta*xi2  );
	      
	    case 3:  // third order
	      mu2.set(6,  eta2     );
	      mu2.set(5,  xi2      );
	      mu2.set(4,  xi*eta   );		    
	      
	    case 2:  // second order		    
	      mu2.set(3, eta       );
	      mu2.set(2, xi        );
	      
	    case 1:  // first order
	      mu2.set(1, 1.0       );
	      
	      break;		    
	    }
	  
	  // Finally, convert monomials to Legendre Polys
	  // on the two adjacent triangle
	  for (int k=1; k<=NumBasisComps; k++)
	    {
	      double tmp1 = 0.0;
	      double tmp2 = 0.0;
	      for (int j=1; j<=k; j++)
		{  
		  tmp1 = tmp1 + Mmat[k-1][j-1]*mu1.get(j);
		  tmp2 = tmp2 + Mmat[k-1][j-1]*mu2.get(j);
		}
	      
	      EdgeData.phi_left->set(i,m,k,  tmp1 );
	      EdgeData.phi_right->set(i,m,k, tmp2 );
	    }
	}
    }
  
}
