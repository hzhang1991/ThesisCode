#include "meshdefs.h"
#include "mesh.h"
#include <math.h>

//
// Post-processing AFTER all boundary and edge information is computed
//
void MeshPostProcess3(const int sub_factor,
		      char*& GridType, 
                      char*& outputdir,
		      mesh& Mesh)
{
  const double xlft = Mesh.get_LowerLeft(1);
  const double ybot = Mesh.get_LowerLeft(2);
  const double xrgt = Mesh.get_UpperRight(1);
  const double ytop = Mesh.get_UpperRight(2);
  const double TOL = 1.0e-12;

  const double oxscale = 1.0/(xrgt-xlft);
  const double oyscale = 1.0/(ytop-ybot);
  const double otmp = 10.0/Mesh.get_shortest_edge();

  const int NumBndNodes = Mesh.get_NumBndNodes();
  int NumPrimBndNodes = 0;
  int NumRdntBndNodes = 0;

  int* prim_tmp = new int[NumBndNodes];
  int* rdnt_tmp = new int[NumBndNodes];

  double* prim_tmp_val = new double[NumBndNodes];
  double* rdnt_tmp_val = new double[NumBndNodes];

  for (int i=1; i<=NumBndNodes; i++)
    {
      const    int k = Mesh.get_bnd_node(i);
      const double x = Mesh.get_node(k,1);
      const double y = Mesh.get_node(k,2);
      
      const bool is_lft = fabs(x-xlft)<TOL;
      const bool is_rgt = fabs(x-xrgt)<TOL;
      const bool is_bot = fabs(y-ybot)<TOL;
      const bool is_top = fabs(y-ytop)<TOL;

      const bool is_corner = ((is_lft && is_bot) ||
			      (is_lft && is_top)) ||
				       ((is_rgt && is_bot) ||
					(is_rgt && is_top));      
      
      const double tmpx = oxscale*(x-xlft);
      const double tmpy = oyscale*(y-ybot);
      double z;

      if (is_corner)
	{ z = 0.0; }
      else
	{
	  if (is_lft || is_rgt)
	    { z = otmp*tmpy; }
	  else if (is_bot || is_top)
	    { z = tmpx; }
	  else
	    {
	      printf(" ERROR in MeshPostProcess3.cpp: \n");
	      printf("  found a boundary subnode that is not on the boundary.\n");
	      printf("\n");
	      exit(1);
	    }
	}

      if (is_lft || is_bot)
	{	  
	  rdnt_tmp[NumRdntBndNodes] = k;	  
	  rdnt_tmp_val[NumRdntBndNodes] = z;
	  NumRdntBndNodes = NumRdntBndNodes + 1;
	}
      else
	{	  
	  prim_tmp[NumPrimBndNodes] = k;
	  prim_tmp_val[NumPrimBndNodes] = z;
	  NumPrimBndNodes = NumPrimBndNodes + 1;
	}
    }

  int* prim_bnd_node = new int[NumPrimBndNodes];
  int* rdnt_bnd_node = new int[NumRdntBndNodes];

  double* prim_bnd_node_val = new double[NumPrimBndNodes];
  double* rdnt_bnd_node_val = new double[NumRdntBndNodes];

  for (int i=0; i<NumPrimBndNodes; i++)
    {
      prim_bnd_node[i] = prim_tmp[i];
      prim_bnd_node_val[i] = prim_tmp_val[i];
    }

  for (int i=0; i<NumRdntBndNodes; i++)
    {
      rdnt_bnd_node[i] = rdnt_tmp[i];
      rdnt_bnd_node_val[i] = rdnt_tmp_val[i];
    }

  void QuickSort(double*& a, int*& index, int lo, int hi);
  QuickSort(prim_bnd_node_val,prim_bnd_node,0,NumPrimBndNodes-1);
  QuickSort(rdnt_bnd_node_val,rdnt_bnd_node,0,NumRdntBndNodes-1);

  // Create node_map
  const int NumPhysNodes = Mesh.get_NumPhysNodes();
  iTensor1 node_map(NumPhysNodes);
  for (int i=1; i<=NumPhysNodes; i++)
    {  
      node_map.set(i, i );  
    }

  {
    int j=0;
    for (int k=0; k<NumRdntBndNodes; k++)
      {
	bool mfound = false;
	while(!mfound)
	  {
	    if (fabs(rdnt_bnd_node_val[k]-prim_bnd_node_val[j])<TOL)
	      {
		mfound = true;		
		node_map.set( rdnt_bnd_node[k], prim_bnd_node[j] );
	      }
	    else
	      {
		j=j+1;
		if (j==NumPrimBndNodes)
		  {
		    printf("\n");
		    printf(" ERROR in MeshPostProcess3.cpp:\n");
		    printf("Could not match primary and redundant boundary node\n");
		    printf("k = %4i,  rdnt_bnd_node[%4i] = %4i, rdnt_bnd_node_val[%4i] = %10.2e\n",
			   k,k,rdnt_bnd_node[k],k,rdnt_bnd_node_val[k]);
		    printf("\n");
		    exit(1);
		  }
	      }
	  }
      }
  }

  char fname1[1024];
  snprintf(fname1,1024,"%s/mesh_node_map.dat",outputdir);
  FILE* write1 = fopen(fname1,"w");
  for (int i=1; i<=NumPhysNodes; i++)
    {  fprintf(write1,"%8i\n",node_map.get(i));  }
  fclose(write1);

  char fname2[1024];
  snprintf(fname2,1024,"%s/mesh_rdnt_node.dat",outputdir);
  FILE* write2 = fopen(fname2,"w");
  fprintf(write2,"%8i\n",NumRdntBndNodes);
  for (int i=0; i<NumRdntBndNodes; i++)
    {  fprintf(write2,"%8i  %8i\n",
	       rdnt_bnd_node[i],
	       node_map.get(rdnt_bnd_node[i]));  }
  fclose(write2);

  // Finished
  delete prim_tmp;
  delete rdnt_tmp;
  delete prim_tmp_val;
  delete rdnt_tmp_val;
  delete prim_bnd_node;
  delete rdnt_bnd_node;
  delete prim_bnd_node_val;
  delete rdnt_bnd_node_val;

  // If necessary, apply procedure to the sub-mesh
  if (sub_factor>1)
    {
      void SubMeshPostProcess3(char*& outputdir,
			       mesh& Mesh);
      SubMeshPostProcess3(outputdir,
			  Mesh);
    }
}

void SubMeshPostProcess3(char*& outputdir,
			 mesh& Mesh)
{
  const double xlft = Mesh.get_LowerLeft(1);
  const double ybot = Mesh.get_LowerLeft(2);
  const double xrgt = Mesh.get_UpperRight(1);
  const double ytop = Mesh.get_UpperRight(2);
  const double TOL = 1.0e-12;

  const double oxscale = 1.0/(xrgt-xlft);
  const double oyscale = 1.0/(ytop-ybot);
  const double otmp = 10.0/Mesh.get_shortest_edge();

  const int SubNumBndNodes = Mesh.get_SubNumBndNodes();
  int SubNumPrimBndNodes = 0;
  int SubNumRdntBndNodes = 0;

  int* prim_tmp = new int[SubNumBndNodes];
  int* rdnt_tmp = new int[SubNumBndNodes];

  double* prim_tmp_val = new double[SubNumBndNodes];
  double* rdnt_tmp_val = new double[SubNumBndNodes];

  for (int i=1; i<=SubNumBndNodes; i++)
    {
      const    int k = Mesh.get_sub_bnd_node(i);
      const double x = Mesh.get_sub_node(k,1);
      const double y = Mesh.get_sub_node(k,2);
      
      const bool is_lft = fabs(x-xlft)<TOL;
      const bool is_rgt = fabs(x-xrgt)<TOL;
      const bool is_bot = fabs(y-ybot)<TOL;
      const bool is_top = fabs(y-ytop)<TOL;

      const bool is_corner = ((is_lft && is_bot) ||
			      (is_lft && is_top)) ||
				       ((is_rgt && is_bot) ||
					(is_rgt && is_top));      
      
      const double tmpx = oxscale*(x-xlft);
      const double tmpy = oyscale*(y-ybot);
      double z;

      if (is_corner)
	{ z = 0.0; }
      else
	{
	  if (is_lft || is_rgt)
	    { z = otmp*tmpy; }
	  else if (is_bot || is_top)
	    { z = tmpx; }
	  else
	    {
	      printf(" ERROR in SubMeshPostProcess3.cpp: \n");
	      printf("  found a boundary subnode that is not on the boundary.\n");
	      printf("\n");
	      exit(1);
	    }
	}

      if (is_lft || is_bot)
	{	  
	  rdnt_tmp[SubNumRdntBndNodes] = k;	  
	  rdnt_tmp_val[SubNumRdntBndNodes] = z;
	  SubNumRdntBndNodes = SubNumRdntBndNodes + 1;
	}
      else
	{	  
	  prim_tmp[SubNumPrimBndNodes] = k;
	  prim_tmp_val[SubNumPrimBndNodes] = z;
	  SubNumPrimBndNodes = SubNumPrimBndNodes + 1;
	}
    }

  int* prim_bnd_node = new int[SubNumPrimBndNodes];
  int* rdnt_bnd_node = new int[SubNumRdntBndNodes];

  double* prim_bnd_node_val = new double[SubNumPrimBndNodes];
  double* rdnt_bnd_node_val = new double[SubNumRdntBndNodes];

  for (int i=0; i<SubNumPrimBndNodes; i++)
    {
      prim_bnd_node[i] = prim_tmp[i];
      prim_bnd_node_val[i] = prim_tmp_val[i];
    }

  for (int i=0; i<SubNumRdntBndNodes; i++)
    {
      rdnt_bnd_node[i] = rdnt_tmp[i];
      rdnt_bnd_node_val[i] = rdnt_tmp_val[i];
    }

  void QuickSort(double*& a, int*& index, int lo, int hi);
  QuickSort(prim_bnd_node_val,prim_bnd_node,0,SubNumPrimBndNodes-1);
  QuickSort(rdnt_bnd_node_val,rdnt_bnd_node,0,SubNumRdntBndNodes-1);

  // Create node_map
  const int SubNumPhysNodes = Mesh.get_SubNumPhysNodes();
  iTensor1 node_map(SubNumPhysNodes);  
  for (int i=1; i<=SubNumPhysNodes; i++)
    {  
      node_map.set(i, i );  
    }

  {
    int j=0;
    for (int k=0; k<SubNumRdntBndNodes; k++)
      {
	bool mfound = false;
	while(!mfound)
	  {
	    if (fabs(rdnt_bnd_node_val[k]-prim_bnd_node_val[j])<TOL)
	      {
		mfound = true;		
		node_map.set( rdnt_bnd_node[k], prim_bnd_node[j] );
	      }
	    else
	      {
		j=j+1;
		if (j==SubNumPrimBndNodes)
		  {
		    printf("\n");
		    printf(" ERROR in SubMeshPostProcess3.cpp:\n");
		    printf("Could not match primary and redundant boundary node\n");
		    printf("k = %4i,  rdnt_bnd_node[%4i] = %4i, rdnt_bnd_node_val[%4i] = %10.2e\n",
			   k,k,rdnt_bnd_node[k],k,rdnt_bnd_node_val[k]);
		    printf("\n");
		    exit(1);
		  }
	      }
	  }
      }
  }

  char fname1[1024];
  snprintf(fname1,1024,"%s/submesh_node_map.dat",outputdir);
  FILE* write1 = fopen(fname1,"w");
  for (int i=1; i<=SubNumPhysNodes; i++)
    {  fprintf(write1,"%8i\n",node_map.get(i));  }
  fclose(write1);

  char fname2[1024];
  snprintf(fname2,1024,"%s/submesh_rdnt_node.dat",outputdir);
  FILE* write2 = fopen(fname2,"w");
  fprintf(write2,"%8i\n",SubNumRdntBndNodes);
  for (int i=0; i<SubNumRdntBndNodes; i++)
    {  fprintf(write2,"%8i  %8i\n",
	       rdnt_bnd_node[i],
	       node_map.get(rdnt_bnd_node[i]));  }
  fclose(write2);

  // Finished
  delete prim_tmp;
  delete rdnt_tmp;
  delete prim_tmp_val;
  delete rdnt_tmp_val;
  delete prim_bnd_node;
  delete rdnt_bnd_node;
  delete prim_bnd_node_val;
  delete rdnt_bnd_node_val;
}
