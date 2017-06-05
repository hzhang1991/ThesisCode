#include "../defs.h"

void SetBC_SubGrid(iTensor2 leftBC, iTensor2 rightBC, double tn, double t, 
		   dTensor2 node, dTensorBC3 aux, dTensorBC4 HermiteCoeffs,
		   dTensorBC3& q_sub)
{
    int i,me,k;
    int mx_sub   = q_sub.getsize(1);
    int meqn_sub = q_sub.getsize(2);
    int kmax_sub = q_sub.getsize(3);
    int mbc_sub  = q_sub.getmbc();
    int mx       = HermiteCoeffs.getsize(1);
    int meqn     = HermiteCoeffs.getsize(2);
    int kmax     = HermiteCoeffs.getsize(3);    
    int mbc      = HermiteCoeffs.getmbc();
    int istart,iend;
    double a0,a1,a2,a3;
    dTensorBC3 qtmp(mx,meqn_sub,kmax,mbc);
    void L2Project(int mopt, int istart, int iend, dTensor2 node,
		   dTensorBC3 qin, dTensorBC3 auxin,  dTensorBC3& Fout,
		   void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&));
    void Map_MainToSub(dTensor1 xpts, dTensor2 q_main, dTensor2 aux_main,
		       dTensor2& q_sub);

    // Left boundary
    if (leftBC.get(mbc_sub,2)==1)
    {
        // Set starting and ending indeces needed by L2Project
	istart = leftBC.get(mbc_sub,1);
	iend   = leftBC.get(1,1);

	// Get solution at the current time
	for (k=1; k<=kmax; k++)
	  for (me=1; me<=meqn; me++)
	    for (i=istart; i<=iend; i++)
	      {
		a0 = HermiteCoeffs.get(i,me,k,1);
		a1 = HermiteCoeffs.get(i,me,k,2);
		a2 = HermiteCoeffs.get(i,me,k,3);
		a3 = HermiteCoeffs.get(i,me,k,4);

		qtmp.set(i,me,k, a0 + a1*(t-tn) + a2*pow(t-tn,2)
			 + a3*pow(t-tn,3) );
	      }

	// Use L2Project with Map_MainToSub and temporarily store
	// in the array qtmp
	L2Project(0,istart,iend,node,qtmp,aux,qtmp,&Map_MainToSub);
	
	// Assign qtmp values to the q_sub array
	for (k=1; k<=kmax_sub; k++)
	  for (me=1; me<=meqn_sub; me++)
	    for (i=1; i<=mbc_sub; i++)
	      {
		q_sub.set(1-i,me,k,  qtmp.get(leftBC.get(i,1),me,k) );
	      }
    }

    // Right boundary
    if (rightBC.get(1,2)==1)
    {
        // Set starting and ending indeces needed by L2Project
	istart = rightBC.get(1,1);
	iend   = rightBC.get(mbc_sub,1);

	// Get solution at the current time
	for (k=1; k<=kmax; k++)
	  for (me=1; me<=meqn; me++)
	    for (i=istart; i<=iend; i++)
	      {
		a0 = HermiteCoeffs.get(i,me,k,1);
		a1 = HermiteCoeffs.get(i,me,k,2);
		a2 = HermiteCoeffs.get(i,me,k,3);
		a3 = HermiteCoeffs.get(i,me,k,4);

		qtmp.set(i,me,k, a0 + a1*(t-tn) + a2*pow(t-tn,2)
			 + a3*pow(t-tn,3) );
	      }
	
	// Use L2Project with Map_MainToSub and temporarily store
	// in the array qtmp
	L2Project(0,istart,iend,node,qtmp,aux,qtmp,&Map_MainToSub);
	
	// Assign qtmp values to the q_sub array
	for (k=1; k<=kmax_sub; k++)
	  for (me=1; me<=meqn_sub; me++)
	    for (i=1; i<=mbc_sub; i++)
	      {
		q_sub.set(mx_sub+i,me,k, qtmp.get(rightBC.get(i,1),me,k) );
	      }
    }

}
