#include "dogdefs.h"
#include "mesh.h"
#include "edge_data_Unst.h"

void BiCGStab(const int MaxIters,
        const double TOL,
        const mesh& Mesh,
        const edge_data_Unst& EdgeData,
        const dTensor2& rhs,
        dTensor2& phi,
        void (*MatMult)(const mesh&,
            const edge_data_Unst&,
            const dTensor2&,
            dTensor2&))
{
    double DotProd(const dTensor2& avec,
            const dTensor2& bvec);
    const int NumPhysElems = phi.getsize(1);
    const int kmax         = phi.getsize(2);
    dTensor2 Aphi(NumPhysElems,kmax);
    dTensor2 res(NumPhysElems,kmax);
    dTensor2 res_star(NumPhysElems,kmax);
    dTensor2 p_res(NumPhysElems,kmax);
    dTensor2 v_res(NumPhysElems,kmax);
    dTensor2 s_res(NumPhysElems,kmax);
    dTensor2 t_res(NumPhysElems,kmax);

    int NumIters  = 0;
    int mflag = 0;

    double rhs_norm = sqrt(DotProd(rhs,rhs));
    if (fabs(rhs_norm)<=1.0e-12)
    { rhs_norm = 1.0; }

    MatMult(Mesh,EdgeData,phi,Aphi);
    for (int i=1; i<=NumPhysElems; i++)
        for (int k=1; k<=kmax; k++)
        {
            res.set(i,k, rhs.get(i,k) - Aphi.get(i,k) ); 
        }

    double rel_res_norm = sqrt(DotProd(res,res))/rhs_norm;
    if(rel_res_norm < TOL)
    {  mflag = 1; }

    for (int i=1; i<=NumPhysElems; i++)
        for (int k=1; k<=kmax; k++)
        {
            res_star.set(i,k, res.get(i,k) );
        }

    double omega = 1;
    double rho;
    double rho1;
    double alpha;
    double beta;

    NumIters = 1;
    while(mflag==0)
    {
        rho = DotProd(res_star,res);

        if (NumIters>1)
        {
            beta  = ( rho/rho1 )*( alpha/omega );
            for (int i=1; i<=NumPhysElems; i++)
                for (int k=1; k<=kmax; k++)
                {
                    p_res.set(i,k, res.get(i,k) 
                            + beta*( p_res.get(i,k) - omega*v_res.get(i,k) ) );
                }
        }
        else
        {
            for (int i=1; i<=NumPhysElems; i++)
                for (int k=1; k<=kmax; k++)
                {
                    p_res.set(i,k, res.get(i,k) );
                }
        }

        MatMult(Mesh,EdgeData,p_res,v_res);
        alpha = rho / DotProd(res_star,v_res);

        for (int i=1; i<=NumPhysElems; i++)
            for (int k=1; k<=kmax; k++)
            {
                s_res.set(i,k, res.get(i,k) - alpha*v_res.get(i,k) );
            }
        double norm_s = sqrt(DotProd(s_res,s_res));

        if (norm_s<TOL)
        {
            for (int i=1; i<=NumPhysElems; i++)
                for (int k=1; k<=kmax; k++)
                {
                    double tmp = phi.get(i,k);
                    phi.set(i,k, tmp + alpha*p_res.get(i,k) );
                }
            rel_res_norm = norm_s/rhs_norm;
            mflag = 1;
        }
        else
        {
            MatMult(Mesh,EdgeData,s_res,t_res);
            omega = DotProd(t_res,s_res)/DotProd(t_res,t_res);

            for (int i=1; i<=NumPhysElems; i++)
                for (int k=1; k<=kmax; k++)
                {
                    double tmp = phi.get(i,k);
                    phi.set(i,k, tmp + alpha*p_res.get(i,k) + omega*s_res.get(i,k) );
                }

            for (int i=1; i<=NumPhysElems; i++)
                for (int k=1; k<=kmax; k++)
                {
                    res.set(i,k, s_res.get(i,k) - omega*t_res.get(i,k) );
                }

            rel_res_norm = sqrt(DotProd(res,res))/rhs_norm;
            if(rel_res_norm < TOL)
            {  mflag = 1;  }

            printf("   NumIters = %i,   rel_res_norm = %e\n",NumIters,rel_res_norm);

            if (NumIters==MaxIters)
            {  mflag = 1;  }

            if (mflag==0)
            {
                rho1 = rho;
                NumIters = NumIters+1;
            }
        }
    }

    printf("  |---------------------------\n");
    printf("  | BiCGStab results:\n");
    printf("  |---------------------------\n");
    printf("  |  MaxIters = %i\n",MaxIters);
    printf("  |       TOL = %e\n",TOL);
    printf("  |  NumIters = %i\n",NumIters);
    printf("  |  residual = %e\n",rel_res_norm);
    printf("  |---------------------------\n");
    printf("\n");
}


double DotProd(const dTensor2& avec,
        const dTensor2& bvec)
{
    const int NumPhysElems = avec.getsize(1);
    const int kmax = avec.getsize(2);
    double prod = 0.0;

    for (int i=1; i<=NumPhysElems; i++)
        for (int k=1; k<=kmax; k++)
        {
            prod = prod + avec.get(i,k)*bvec.get(i,k);
        }

    return prod;
}
