#include "../defs.h"

void GetHermiteCoeffs(int method[],int istart, int iend, double dt, double tn,
		      dTensorBC3 qn, dTensorBC3 auxn, dTensorBC3 qn1, 
		      dTensorBC3 auxn1, dTensor2 node, dTensorBC4& coeffs)
{
    int mx   = qn.getsize(1);
    int meqn = qn.getsize(2);
    int kmax = qn.getsize(3);
    int maux = auxn.getsize(2);
    int mbc  = qn.getmbc();
    dTensorBC3  Ln(mx,meqn,kmax,mbc);
    dTensorBC3 Ln1(mx,meqn,kmax,mbc);
    int i,m,k;
    double Qn_tmp,Fn_tmp,Qn1_tmp,Fn1_tmp;
    void ConstructL_HMM_short(int,int,int[],dTensor2,dTensorBC3,
			      dTensorBC3,dTensorBC3&,
			      void (*FluxFunc)(dTensor1,dTensor2,dTensor2,dTensor2&),
			      void (*SetWaveSpd)(dTensor1,dTensor2,dTensor2,dTensor2,
						 dTensor2,double&,double&),
			      void (*SourceTermFunc)(dTensor1,dTensor2,
						     dTensor2,dTensor2&));
    void FluxFunc(dTensor1,dTensor2,dTensor2,dTensor2&);
    void SetWaveSpd(dTensor1,dTensor2,dTensor2,dTensor2,dTensor2,double&,double&);
    void SourceTermFunc(dTensor1,dTensor2,dTensor2,dTensor2&);

    // Get q_t = L(q) -- time tn
    ConstructL_HMM_short(istart,iend,method,node,auxn,qn,Ln,&FluxFunc,
			 &SetWaveSpd,&SourceTermFunc);

    // Get q_t = L(q) -- time tn1
    ConstructL_HMM_short(istart,iend,method,node,auxn1,qn1,Ln1,&FluxFunc,
			 &SetWaveSpd,&SourceTermFunc);

    // Construct coefficients for the Hermite Interpolating Polynomial
    for (k=1; k<=kmax; k++)
      for (m=1; m<=meqn; m++)
	for (i=istart; i<=iend; i++)
	  {
	    Qn_tmp  = qn.get(i,m,k);
	    Fn_tmp  = Ln.get(i,m,k);
	    Qn1_tmp = qn1.get(i,m,k);
	    Fn1_tmp = Ln1.get(i,m,k);

	    coeffs.set(i,m,k,1,  Qn_tmp );
	    coeffs.set(i,m,k,2,  Fn_tmp );
	    coeffs.set(i,m,k,3, (3.0*(Qn1_tmp-Qn_tmp) - dt*(2.0*Fn_tmp+Fn1_tmp))/(pow(dt,2)) );
	    coeffs.set(i,m,k,4, (2.0*(Qn_tmp-Qn1_tmp) + dt*(Fn_tmp+Fn1_tmp))/(pow(dt,3)) );
	  }
}
