#include "../defs.h"

// Advance forward one stage (Runge-Kutta method)
void AdvanceOneStage_RK(double alpha1, double alpha2, double beta, double dt,
			int method[], dTensor2 node, dTensorBC3 qn,
			dTensorBC3 aux_old, dTensorBC3 q_old,
			dTensorBC3& aux_new, dTensorBC3& q_new,
			dTensorBC1& smax, iTensor1 map,
			iTensor2 leftBC, iTensor2 rightBC,int mainsub)
{
    int mx     = q_old.getsize(1);
    int meqn   = q_old.getsize(2);
    int mbc    = q_old.getmbc();
    int maux   = aux_old.getsize(2);
    dTensorBC3   Lrhs(mx,meqn,method[1],mbc);    
    void ConstructL_HMM(int[],dTensor2,dTensorBC3,dTensorBC3,
			dTensorBC3&,dTensorBC1&,iTensor1,
			iTensor2,iTensor2,int,
			void (*FluxFunc)(dTensor1,dTensor2,dTensor2,dTensor2&),
			void (*SetWaveSpd)(dTensor1,dTensor2,dTensor2,dTensor2,
					   dTensor2,double&,double&),
			void (*SourceTermFunc)(dTensor1,dTensor2,dTensor2,dTensor2&));
    void ApplyLimiter(dTensor2,dTensorBC3,dTensorBC3&,
		      void (*ProjectRightEig)(dTensor1,dTensor1,dTensor2,dTensor2&),
		      void (*ProjectLeftEig)(dTensor1,dTensor1,dTensor2,dTensor2&));    
    void UpdateSoln_HMM(double,double,double,double,dTensor2,
			dTensorBC3,dTensorBC3,dTensorBC3,dTensorBC3&);

    // MainGrid
    void BeforeStep(dTensor2,dTensorBC3&,dTensorBC3&);
    void AfterStep(dTensor2,dTensorBC3&,dTensorBC3&);
    void FluxFunc(dTensor1,dTensor2,dTensor2,dTensor2&);
    void SetWaveSpd(dTensor1,dTensor2,dTensor2,dTensor2,dTensor2,double&,double&);
    void SourceTermFunc(dTensor1,dTensor2,dTensor2,dTensor2&);
    void ProjectRightEig(dTensor1,dTensor1,dTensor2,dTensor2&);
    void ProjectLeftEig(dTensor1,dTensor1,dTensor2,dTensor2&);

    // SubGrid
    void BeforeStep_sub(dTensor2,dTensorBC3&,dTensorBC3&);
    void AfterStep_sub(dTensor2,dTensorBC3&,dTensorBC3&);
    void FluxFunc_sub(dTensor1,dTensor2,dTensor2,dTensor2&);
    void SetWaveSpd_sub(dTensor1,dTensor2,dTensor2,dTensor2,dTensor2,double&,double&);
    void SourceTermFunc_sub(dTensor1,dTensor2,dTensor2,dTensor2&);
    void ProjectRightEig_sub(dTensor1,dTensor1,dTensor2,dTensor2&);
    void ProjectLeftEig_sub(dTensor1,dTensor1,dTensor2,dTensor2&);
    

    // --------------------------------------------------------------
    // One stage of RK: advances (aux_old, q_old) to (aux_new, q_new)
    // --------------------------------------------------------------

    if (mainsub==1) // MainGrid
    {
	// Set information necessary before stage starts
	BeforeStep(node,aux_old,q_old);
	
	// Method of lines (q_t = L(q)) -- construct RHS L(q)
	ConstructL_HMM(method,node,aux_old,q_old,Lrhs,smax,
		       map,leftBC,rightBC,mainsub,&FluxFunc,&SetWaveSpd,
		       &SourceTermFunc);

	// Update solution based on above computed L(q)
	UpdateSoln_HMM(alpha1,alpha2,beta,dt,node,qn,q_old,Lrhs,q_new);
	
	// Apply limiters
	if (method[3]==1 && method[1]>1)
	  {  ApplyLimiter(node,aux_old,q_new,&ProjectRightEig,&ProjectLeftEig);  }

	// Set information necessary after stage ends
	AfterStep(node,aux_new,q_new);
    }
    else // SubGrid
    {
	// Set information necessary before stage starts
	BeforeStep_sub(node,aux_old,q_old);
	
	// Method of lines (q_t = L(q)) -- construct RHS L(q)
	ConstructL_HMM(method,node,aux_old,q_old,Lrhs,smax,
		       map,leftBC,rightBC,mainsub,&FluxFunc_sub,&SetWaveSpd_sub,
		       &SourceTermFunc_sub);

	// Update solution based on above computed L(q)
	UpdateSoln_HMM(alpha1,alpha2,beta,dt,node,qn,q_old,Lrhs,q_new);
	
	// Apply limiters
	if (method[3]==1 && method[1]>1)
	  {  ApplyLimiter(node,aux_old,q_new,&ProjectRightEig_sub,&ProjectLeftEig_sub);  }

	// Set information necessary after stage ends
	AfterStep_sub(node,aux_new,q_new);
    }

}
