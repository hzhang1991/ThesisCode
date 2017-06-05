#include "../defs.h"

void SetValues_MainGrid(iTensor1 map, dTensor2 node_sub, dTensorBC3 aux_sub,
			dTensorBC3 q_sub, dTensor2 node, dTensorBC3 aux, 
			dTensorBC3& q)
{
    int i,me,k;
    int mx_sub   = q_sub.getsize(1);
    int meqn_sub = q_sub.getsize(2);
    int kmax_sub = q_sub.getsize(3);
    int mx       = q.getsize(1);
    int meqn     = q.getsize(2);
    int kmax     = q.getsize(3);
    int mbc_sub  = q_sub.getmbc();
    int mbc      = q.getmbc();
    int maux_sub = aux_sub.getsize(2);
    int istart,iend;
    dTensorBC3 qtmp(mx_sub,meqn,kmax_sub,mbc_sub);
    void L2Project(int mopt, int istart, int iend, dTensor2 node,
		   dTensorBC3 qin, dTensorBC3 auxin,  dTensorBC3& Fout,
		   void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&));
    void Map_SubToMain(dTensor1 xpts, dTensor2 q_sub, dTensor2 aux_sub,
		       dTensor2& q_main);

    // Set starting and ending indeces needed by L2Project
    istart = 1;
    iend   = mx_sub;

    // Use L2Project with Map_MainToSub and temporarily store
    // in the array qtmp
    L2Project(0,istart,iend,node_sub,q_sub,aux_sub,qtmp,&Map_SubToMain);
	
    // Assign qtmp values to the q array
    for (k=1; k<=kmax; k++)
      for (me=1; me<=meqn; me++)
	for (i=1; i<=mx_sub; i++)
	  {
	    q.set(map.get(i),me,k,  qtmp.get(i,me,k) );
	  }

}
