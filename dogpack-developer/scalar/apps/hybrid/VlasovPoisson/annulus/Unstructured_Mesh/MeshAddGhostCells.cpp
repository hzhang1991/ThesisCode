#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include "meshdefs.h"

inline double modxy(const double& xy,
		    const double& xylow,
		    const double& xyhigh,
		    const double& lxy)
{
  double xyval;

  if (xy>xyhigh)
    { xyval = xy-lxy; }
  else if (xy<xylow)
    { xyval = xy+lxy; }
  else
    { xyval = xy; }

  return xyval;
}
  
inline int find_center(const int& num_centers,
		       const int* centers_index,
		       const double* centers_xy,
		       const double& cent)
{
  bool found = 0;
  int k = -1;
  int val;

  while(!found)
    {
      k=k+1;
      if (k==num_centers)
	{
	  printf("\n");
	  printf(" Error in find_center: did not find center ...\n");
	  printf("   cent = %e\n",cent);
	  printf("\n");
	  exit(1);
	}

      if (fabs(cent-centers_xy[k])<1.0e-12)
	{
	  val = centers_index[k];
	  found = 1;
	}
    }

  return val;
}

//
// Add ghost cells
//
void MeshAddGhostCells(const string& GridType, cart2d& cart2dinfo,
		       int& numpts, int& numtri, int& numghost, int& num_ext_node,
		       point*& p, triangle*& t, double*& area, double*& cdual, 
		       int*& ghost_link, int*& proper_ghostcell,
		       int*& ext_node_link,
		       double (*SignedDistance)(point))
{
  int i,j;
  int k = 0;
  int mfound = 0;
  int jstore[4];
  double x,y;
  int tmp[4];
  point ptmp;
  point A,B,C,nvec,Anew,v12,v23,v31;
  double BCnorm,nvec_norm;
  int tmps,temp,kcheck;
  point alpha,beta;
  double alpha_norm,beta_norm,theta_alpha,theta_beta;
  int jopp;
  
  point pnew[2*numpts];
  triangle tnew[2*numtri];
  double area_new[2*numtri];
  double cdual_new[2*numpts];
  int ghost_link_new[2*numtri];
  int numghost_new = 0;
  int numpts_new = numpts;
  int numtri_new = numtri;
  int ext_new[2*numtri];
  int num_ext = 0;

  // for periodic boundary conditions
  int* centers_index_tmp = new int[2*numtri];
  double* centers_xy_tmp = new double[2*numtri];
  int num_centers = 0;
  const double olx = 1.0/(cart2dinfo.xhigh-cart2dinfo.xlow);
  const double oly = 1.0/(cart2dinfo.yhigh-cart2dinfo.ylow);

  // Compute the minimal distance parameter hmin
  double min_parea = 1.0e10;
  for (i=0; i<numtri; i++)
    {
      if (min_parea > area[i])
	{ min_parea = area[i]; }
    }
  double hmin = sqrt(2.0*min_parea)/5.0e0;

  // Store old t,area,p,cdual in new versions of each
  for (i=0; i<numtri; i++)
    {
      tnew[i].n1 = t[i].n1;
      tnew[i].n2 = t[i].n2;
      tnew[i].n3 = t[i].n3;
      
      area_new[i] = area[i];
    }

  for (i=0; i<numpts; i++)
    {
      pnew[i].x = p[i].x;
      pnew[i].y = p[i].y;

      cdual_new[i] = cdual[i];
    }

  // Loop over every element in order to determine
  // if a ghost cell should be created
  for (i=0; i<numtri; i++)
    {
      mfound = 0;
      jstore[1] = 0;
      jstore[2] = 0;
      jstore[3] = 0;

      // Look at ith element
      tmp[1] = t[i].n1;
      tmp[2] = t[i].n2;
      tmp[3] = t[i].n3;

      // Check how many nodes of ith element lie on the boundary (either 0, 1, 2, or 3)
      // If there are 0 or 1 nodes on the bnd, then do nothing and proceed to next element
      for (j=1; j<=3; j++)
	{
	  ptmp.x = p[tmp[j]].x;
	  ptmp.y = p[tmp[j]].y;
	  
	  if (fabs(SignedDistance(ptmp))<=hmin)
	    {
	      mfound = mfound + 1;
	      jstore[mfound] = j;
	    }
	  else
	    {
	      A.x = ptmp.x;
	      A.y = ptmp.y;
	      jopp = tmp[j];
	    }
	}    

      // Exactly three nodes of a triangle are on the boundary
      // NOTE: 3 nodes on boundary <==> exactly two edges are boundary edges
      //   (assuming the mesh has more than 1 element and these elements are connected)
      if (mfound==3)
	{
	  kcheck = 0;

	  // check if edge 12 is on the boundary
	  ptmp.x = 0.5*( p[tmp[1]].x + p[tmp[2]].x );
	  ptmp.y = 0.5*( p[tmp[1]].y + p[tmp[2]].y );

	  if (fabs(SignedDistance(ptmp))<=hmin)
	    {
	      kcheck=kcheck+1;
	      
	      jstore[1] = 1;
	      jstore[2] = 2;
	      jstore[3] = 3;

	      A.x = p[tmp[jstore[3]]].x;
	      A.y = p[tmp[jstore[3]]].y;
	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;
	      
	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;

	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
  
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;

	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = tmp[jstore[3]];
	      
	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;

	      // Center information - to be used for periodic BCs
	      num_centers = num_centers + 1;
	      centers_index_tmp[num_centers-1] = i;
	      double x1 = A.x;
	      double x2 = B.x;
	      double x3 = C.x;
	      double y1 = A.y;
	      double y2 = B.y;
	      double y3 = C.y;
	      centers_xy_tmp[num_centers-1] = olx*(onethird*(x1+x2+x3)-cart2dinfo.xlow)
		+ 1000.0*oly*(onethird*(y1+y2+y3)-cart2dinfo.ylow);	      
	      
	      // Compute signed element area, change orientation if necessary
	      v12.x = pnew[tnew[numtri_new-1].n2].x - pnew[tnew[numtri_new-1].n1].x;
	      v12.y = pnew[tnew[numtri_new-1].n2].y - pnew[tnew[numtri_new-1].n1].y;
	      v23.x = pnew[tnew[numtri_new-1].n3].x - pnew[tnew[numtri_new-1].n2].x;
	      v23.y = pnew[tnew[numtri_new-1].n3].y - pnew[tnew[numtri_new-1].n2].y;
	      v31.x = pnew[tnew[numtri_new-1].n1].x - pnew[tnew[numtri_new-1].n3].x;
	      v31.y = pnew[tnew[numtri_new-1].n1].y - pnew[tnew[numtri_new-1].n3].y;
	      
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));

	      if (area_new[numtri_new-1] < 0.0)
		{
		  temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
		}

	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
	    }

	  // check if edge 13 is on the boundary
	  ptmp.x = 0.5*( p[tmp[1]].x + p[tmp[3]].x );
	  ptmp.y = 0.5*( p[tmp[1]].y + p[tmp[3]].y );

	  if (fabs(SignedDistance(ptmp))<=hmin)
	    {
	      kcheck=kcheck+1;
	      
	      jstore[1] = 3;
	      jstore[2] = 1;
	      jstore[3] = 2;

	      A.x = p[tmp[jstore[3]]].x;
	      A.y = p[tmp[jstore[3]]].y;
	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;
	      
	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;

	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
  
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;

	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = tmp[jstore[3]];
	      
	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;

	      // Center information - to be used for periodic BCs
	      num_centers = num_centers + 1;
	      centers_index_tmp[num_centers-1] = i;
	      double x1 = A.x;
	      double x2 = B.x;
	      double x3 = C.x;
	      double y1 = A.y;
	      double y2 = B.y;
	      double y3 = C.y;
	      centers_xy_tmp[num_centers-1] = olx*(onethird*(x1+x2+x3)-cart2dinfo.xlow)
		+ 1000.0*oly*(onethird*(y1+y2+y3)-cart2dinfo.ylow);
	      
	      // Compute signed element area, change orientation if necessary
	      v12.x = pnew[tnew[numtri_new-1].n2].x - pnew[tnew[numtri_new-1].n1].x;
	      v12.y = pnew[tnew[numtri_new-1].n2].y - pnew[tnew[numtri_new-1].n1].y;
	      v23.x = pnew[tnew[numtri_new-1].n3].x - pnew[tnew[numtri_new-1].n2].x;
	      v23.y = pnew[tnew[numtri_new-1].n3].y - pnew[tnew[numtri_new-1].n2].y;
	      v31.x = pnew[tnew[numtri_new-1].n1].x - pnew[tnew[numtri_new-1].n3].x;
	      v31.y = pnew[tnew[numtri_new-1].n1].y - pnew[tnew[numtri_new-1].n3].y;
	      
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));

	      if (area_new[numtri_new-1] < 0.0)
		{
		  temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
		}

	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
	    }

	  // check if edge 23 is on the boundary
	  ptmp.x = 0.5*( p[tmp[2]].x + p[tmp[3]].x );
	  ptmp.y = 0.5*( p[tmp[2]].y + p[tmp[3]].y );

	  if (fabs(SignedDistance(ptmp))<=hmin)
	    {
	      kcheck=kcheck+1;

	      jstore[1] = 2;
	      jstore[2] = 3;
	      jstore[3] = 1;

	      A.x = p[tmp[jstore[3]]].x;
	      A.y = p[tmp[jstore[3]]].y;
	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;
	      
	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;

	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
  
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;

	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = tmp[jstore[3]];
	      
	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;

	      // Center information - to be used for periodic BCs
	      num_centers = num_centers + 1;
	      centers_index_tmp[num_centers-1] = i;
	      double x1 = A.x;
	      double x2 = B.x;
	      double x3 = C.x;
	      double y1 = A.y;
	      double y2 = B.y;
	      double y3 = C.y;
	      centers_xy_tmp[num_centers-1] = olx*(onethird*(x1+x2+x3)-cart2dinfo.xlow)
		+ 1000.0*oly*(onethird*(y1+y2+y3)-cart2dinfo.ylow);
	      
	      // Compute signed element area, change orientation if necessary
	      v12.x = pnew[tnew[numtri_new-1].n2].x - pnew[tnew[numtri_new-1].n1].x;
	      v12.y = pnew[tnew[numtri_new-1].n2].y - pnew[tnew[numtri_new-1].n1].y;
	      v23.x = pnew[tnew[numtri_new-1].n3].x - pnew[tnew[numtri_new-1].n2].x;
	      v23.y = pnew[tnew[numtri_new-1].n3].y - pnew[tnew[numtri_new-1].n2].y;
	      v31.x = pnew[tnew[numtri_new-1].n1].x - pnew[tnew[numtri_new-1].n3].x;
	      v31.y = pnew[tnew[numtri_new-1].n1].y - pnew[tnew[numtri_new-1].n3].y;
	      
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));

	      if (area_new[numtri_new-1] < 0.0)
		{
		  temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
		}

	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
	    }

	  if (kcheck!=2)
	    {
	      cout << endl;
	      cout << " ERROR in MeshAddGhostCells.cpp: " << endl; 
	      cout << "       Found a triangle with 3 nodes on the boundary." << endl;
	      cout << "       This triangle should have exactly 2 edges on the boundary," << endl;
	      cout << "         but in fact it has " << kcheck << endl;
	      cout << endl;
	      exit(1);
	    }
	}
      // Exactly two nodes of a triangle are on the boundary
      // NOTE: 2 nodes on boundary  <==> either exactly one edge is a boundary edge  -or-
      //                                 near a corner and exactly zero edges are a bnd edge
      else if (mfound==2)
	{

	  // Compute edge midpoint
	  ptmp.x = 0.5*( p[tmp[jstore[1]]].x + p[tmp[jstore[2]]].x );
	  ptmp.y = 0.5*( p[tmp[jstore[1]]].y + p[tmp[jstore[2]]].y );

	  // Proceed if edge midpoint is also on the bnd <==> exactly one edge is a bnd edge
	  if (fabs(SignedDistance(ptmp))<=hmin)
	    {
	      if (jstore[1]>jstore[2])
		{
		  tmps = jstore[2];
		  jstore[2] = jstore[1];
		  jstore[1] = tmps;
		}
	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;
	      
	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;

	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
  
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;
	      
	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = jopp;// tmp[jstore[3]];

	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;

	      // Center information - to be used for periodic BCs
	      num_centers = num_centers + 1;
	      centers_index_tmp[num_centers-1] = i;
	      double x1 = A.x;
	      double x2 = B.x;
	      double x3 = C.x;
	      double y1 = A.y;
	      double y2 = B.y;
	      double y3 = C.y;
	      centers_xy_tmp[num_centers-1] = olx*(onethird*(x1+x2+x3)-cart2dinfo.xlow)
		+ 1000.0*oly*(onethird*(y1+y2+y3)-cart2dinfo.ylow);
	      
	      // Compute signed element area, change orientation if necessary
	      v12.x = pnew[tnew[numtri_new-1].n2].x - pnew[tnew[numtri_new-1].n1].x;
	      v12.y = pnew[tnew[numtri_new-1].n2].y - pnew[tnew[numtri_new-1].n1].y;
	      v23.x = pnew[tnew[numtri_new-1].n3].x - pnew[tnew[numtri_new-1].n2].x;
	      v23.y = pnew[tnew[numtri_new-1].n3].y - pnew[tnew[numtri_new-1].n2].y;
	      v31.x = pnew[tnew[numtri_new-1].n1].x - pnew[tnew[numtri_new-1].n3].x;
	      v31.y = pnew[tnew[numtri_new-1].n1].y - pnew[tnew[numtri_new-1].n3].y;
	      
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));

	      if (area_new[numtri_new-1] < 0.0)
		{
		  temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
		}

	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
	    }
	}
    }

  // Resize all appropriate vectors
  delete[] p;
  delete[] t;
  delete[] area;
  delete[] cdual;

  numpts = numpts_new;
  numtri = numtri_new;
  numghost = numghost_new;
  num_ext_node = num_ext;

  p = new point[numpts];
  t = new triangle[numtri];
  area = new double[numtri];
  cdual = new double[numpts];
  ghost_link = new int[numghost];
  proper_ghostcell = new int[numghost];
  ext_node_link = new int[num_ext];

  for (i=0; i<numtri; i++)
    {
      t[i].n1 = tnew[i].n1;
      t[i].n2 = tnew[i].n2;
      t[i].n3 = tnew[i].n3;

      area[i] = area_new[i];
    }

  for (i=0; i<numpts; i++)
    {
      p[i].x = pnew[i].x;
      p[i].y = pnew[i].y;

      cdual[i] = cdual_new[i];
    }

  for (i=0; i<numghost; i++)
    {
      ghost_link[i] = ghost_link_new[i];
    }

  for (i=0; i<numghost; i++)
    {
      proper_ghostcell[i] = 1;
    }

  for (i=0; i<num_ext; i++)
    {
      ext_node_link[i] = ext_new[i];
    }

  // Apply periodic boundary conditions
  int*  centers_index = new int[num_centers];
  double*  centers_xy = new double[num_centers];
  bool* centers_found = new bool[num_centers];
 
  for (i=0; i<num_centers; i++)
    {
      centers_index[i] = centers_index_tmp[i];
      centers_xy[i]    = centers_xy_tmp[i];
      centers_found[i] = 0;
    }
  delete[] centers_index_tmp;
  delete[] centers_xy_tmp;

  const double lx = (cart2dinfo.xhigh-cart2dinfo.xlow);
  const double ly = (cart2dinfo.yhigh-cart2dinfo.ylow);

  for (i=0; i<numghost; i++)
    {
      int n1 = t[numtri-numghost+i].n1;
      int n2 = t[numtri-numghost+i].n2;
      int n3 = t[numtri-numghost+i].n3;
      double x1 = p[n1].x;
      double y1 = p[n1].y;
      double x2 = p[n2].x;
      double y2 = p[n2].y;
      double x3 = p[n3].x;
      double y3 = p[n3].y;
      double xc = modxy(onethird*(x1+x2+x3),cart2dinfo.xlow,cart2dinfo.xhigh,lx);
      double yc = modxy(onethird*(y1+y2+y3),cart2dinfo.ylow,cart2dinfo.yhigh,ly);
      double cent = olx*(xc-cart2dinfo.xlow) + 1000.0*oly*(yc-cart2dinfo.ylow);

      ghost_link[i] = find_center(num_centers,centers_index,centers_xy,cent);
    }

}
