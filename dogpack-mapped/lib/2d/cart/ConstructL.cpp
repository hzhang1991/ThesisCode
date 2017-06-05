//#undef CHECK_BOUNDS // this file is okay, so omit tensor bounds check
#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "L2Project.h"
#include "IntegrateBasis.h"
#include "RiemannSolve.h"
#include "Legendre2d.h"
#include "float.h" // for debugging


// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t = RHS =  -[ f(q,x,y,t)_x + g(q,x,y,t)_y ] + Psi(q,x,y,t)
//


void mapc2p(double& xc,double& yc);
vector<double> returnleft(double x1l,double y1l,double x2l,double y2l,double x3l,double y3l,double x4l,double y4l,int);
vector<double> returnright(double x1l,double y1l,double x2l,double y2l,double x3l,double y3l,double x4l,double y4l,int);



void ConstructL(const dTensorBC4& aux, const dTensorBC4& q,
        dTensorBC4& Lstar, dTensorBC3& smax)
{

    const edge_data& edgeData = Legendre2d::get_edgeData();
    const int space_order = dogParams.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    dTensorBC4 Fm(mx,my,meqn,space_order,mbc);
    dTensorBC4 Fp(mx,my,meqn,space_order,mbc);
    dTensorBC4 Gm(mx,my,meqn,space_order,mbc);
    dTensorBC4 Gp(mx,my,meqn,space_order,mbc);

    // If you need access elsewhere to one of these declarations
    // move it into DogSolverCart2.h rather than copying the
    // declaration. (By using a single declaration we can append
    // arguments with default values without modifying existing
    // calls.) -eaj
    //
    void FluxFunc(const dTensor2& xpts,
            const dTensor2& Q,
            const dTensor2& Aux,
            dTensor3& flux);
    void LstarExtra(const dTensorBC4*,
            const dTensorBC4*,
            dTensorBC4*);
    void ArtificialViscosity(const dTensorBC4* aux, 
            const dTensorBC4* q, 
            dTensorBC4* Lstar);

    // Grid information
    const double xlower = dogParamsCart2.get_xlow();
    const double ylower = dogParamsCart2.get_ylow();
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double xlow=xlower;
    const double ylow=ylower;
   

    // --------------------------------------------------------------------- //
    // Boundary data:
    //    User is now responsible for setting this before calling ConstructL.
    // --------------------------------------------------------------------- //
//  void SetBndValues(dTensorBC4& q, dTensorBC4& aux);
//  SetBndValues(q, aux);
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part I: compute source term
    // --------------------------------------------------------- 
    if ( dogParams.get_source_term()>0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals,
                const dTensor2& auxvals, dTensor2& source);
        L2Project(1-mbc,mx+mbc,1-mbc,my+mbc,
                space_order,space_order,space_order,space_order,
                &q,&aux,&Lstar,&SourceTermFunc);
    }
    else
    {
        Lstar.setall(0.);
    }

    if ( !dogParams.get_flux_term() )
    {  return;  }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part II: compute inter-element interaction fluxes
    // ---------------------------------------------------------
    // 1-direction: loop over interior edges and solve Riemann problems
    /*dTensor1 nvec(2);
    nvec.set(1, 1.0e0 );
    nvec.set(2, 0.0e0 );

#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc); i++)
    {
        dTensor1 Ql(meqn),Qr(meqn);
        dTensor1 Auxl(maux),Auxr(maux);
        dTensor1 Fl(meqn),Fr(meqn);
        dTensor1 xedge(2);

        RiemannSolver riemannSolver(meqn,maux);
        for (int j=(2-mbc); j<=(my+mbc-1); j++)
        {


    double xil1 = xlower + (double(i)-1.0)*dx;
        double yil1 = ylower + (double(j)-1.0)*dy;
        double xil2 = xlower + (double(i))*dx;
        double yil2 = ylower + (double(j-1))*dy; 
        double xil3 = xlower + (double(i))*dx;
        double yil3 = ylower + (double(j))*dy;
        double xil4 = xlower + (double(i)-1)*dx;
        double yil4 = ylower + (double(j))*dy;
        vector<double> x1=returnright(xil1,yil1,xil2,yil2,xil3,yil3,xil4,yil4,2);
        double ax1l = x1.at(2);
        double ay1l = x1.at(3); 
        double ax2l = x1.at(4);
        double ay2l = x1.at(5);
        double ax3l = x1.at(6);
        double ay3l = x1.at(7);


               double x1p = xlow + (double(i+1)-1.0)*dx;
               double y1p = ylow + (double(j+1)-1.0)*dy;
               double x2p = xlow + (double(i+1))*dx;
               double y2p = ylow + (double(j+1-1))*dy; 
               double x3p = xlow + (double(i+1))*dx;
               double y3p = ylow + (double(j+1))*dy;
               double x4p = xlow + (double(i+1)-1)*dx;
               double y4p = ylow + (double(j+1))*dy;
               vector<double> x=returnright(x1p,y1p,x2p,y2p,x3p,y3p,x4p,y4p,2);
               double x1l = x.at(2);
               double y1l = x.at(3); 
               double x2l = x.at(4);
               double y2l = x.at(5);
               double x3l = x.at(6);
               double y3l = x.at(7);

               mapc2p(x1p,y1p);mapc2p(x2p,y2p);mapc2p(x3p,y3p);mapc2p(x4p,y4p);

                nvec.set(1, -(y1p-y4p) );
                nvec.set(2, -(x4p-x1p) ); 
                double mag=0.5/sqrt(nvec.get(1)*nvec.get(1)+nvec.get(2)*nvec.get(2));
                nvec.set(1, -(y1p-y4p)*2.0*mag );
                nvec.set(2, -(x4p-x1p)*2.0*mag ); 

            // ell indexes Riemann point along the edge
            for (int ell=1; ell<=space_order; ell++)
            {
                // Riemann data - q
                for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, 0.0 );
                    Qr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Ql.fetch(m) += edgeData.phi_xl->get(ell,k)*q.get(i-1,j,m,k);
                        Qr.fetch(m) += edgeData.phi_xr->get(ell,k)*q.get(i,j,m,k);
                    }
                }

                // Riemann data - aux
                for (int m=1; m<=maux; m++)
                {
                    Auxl.set(m, 0.0 );
                    Auxr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Auxl.fetch(m) += edgeData.phi_xl->get(ell,k)*aux.get(i-1,j,m,k);
                        Auxr.fetch(m) += edgeData.phi_xr->get(ell,k)*aux.get(i,j,m,k);
                    }
                }

                // Solve Riemann problem
               double s = edgeData.xpts1d->get(ell);
              //  xedge.set(1, xlower + (double(i)-1.0)*dx );
              //  xedge.set(2, ylower + (double(j)-0.5)*dy );
                xedge.set(1, x4p + 0.5*(s+1.0)*(x1p-x4p) );
                xedge.set(2, y4p + 0.5*(s+1.0)*(y1p-y4p) );
               // cout<<"xdiff "<<x1p + 0.5*(1.0)*(x2p-x1p)<<" "<<xlower + (double(i)-0.5)*dx<<endl;
               // cout<<"ydiff "<<y1p + 0.5*(1.0)*(y2p-y1p) <<" "<<ylower + (double(j)-1.0)*dy<<endl;

                const double smax_edge
                    = riemannSolver.solve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);

                // TODO: See comments in all ConstructL methods.  This *should* be a 
                // bug when OpenMP loops are turned on. (-DS)
                //
                // If we reverse the role of i and j in the above for loop,
                // then it will be OK, but the tensors don't run as fast when
                // that's the case.
                //

                smax.set(i-1, j, 1, Max(0.5/mag*smax_edge, smax.get(i-1, j, 1) ) );
                smax.set(i,   j, 1, Max(0.5/mag*smax_edge, smax.get(i,   j, 1) ) );

                //smax.set(i-1, j,1,  Max(abs(smax_edge), smax.get(i, j,1) ) );
                //smax.set(i,   j,1,  Max(abs(smax_edge), smax.get(i+1,   j,1) ) );
                
                // Construct fluxes
                for (int m=1; m<=meqn; m++)
                {
                    Fm.set(i  ,j,m,ell, Fr.get(m)*mag );
                    Fp.set(i-1,j,m,ell, Fl.get(m)*mag );
                }
            }
        }
    }


    // 2-direction: loop over interior edges and solve Riemann problems
    nvec.set(1, 0.0e0 );
    nvec.set(2, 1.0e0 );

#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    {
        dTensor1 Ql(meqn),Qr(meqn);
        dTensor1 Auxl(maux),Auxr(maux);
        dTensor1 xedge(2);
        dTensor1 Fl(meqn),Fr(meqn);

        RiemannSolver riemannSolver(meqn,maux);
        for (int j=(2-mbc); j<=(my+mbc); j++)
            {
               double x1p = xlow + (double(i+1)-1.0)*dx;
               double y1p = ylow + (double(j+1)-1.0)*dy;
               double x2p = xlow + (double(i+1))*dx;
               double y2p = ylow + (double(j+1-1))*dy; 
               double x3p = xlow + (double(i+1))*dx;
               double y3p = ylow + (double(j+1))*dy;
               double x4p = xlow + (double(i+1)-1)*dx;
               double y4p = ylow + (double(j+1))*dy;
               vector<double> x=returnright(x1p,y1p,x2p,y2p,x3p,y3p,x4p,y4p,2);
               double x1l = x.at(2);
               double y1l = x.at(3); 
               double x2l = x.at(4);
               double y2l = x.at(5);
               double x3l = x.at(6);
               double y3l = x.at(7);

               mapc2p(x1p,y1p);mapc2p(x2p,y2p);mapc2p(x3p,y3p);mapc2p(x4p,y4p);
         
               
               nvec.set(1, -(y2p-y1p) );
               nvec.set(2, -(x1p-x2p) );
               double mag=0.5/sqrt(nvec.get(1)*nvec.get(1)+nvec.get(2)*nvec.get(2));
               nvec.set(1, -(y2p-y1p)*2.0*mag );
               nvec.set(2, -(x1p-x2p)*2.0*mag ); 
               

            for (int ell=1; ell<=space_order; ell++)
            {
                // Riemann data - q
                for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, 0.0 );
                    Qr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Ql.fetch(m) += edgeData.phi_yl->get(ell,k)*q.get(i,j-1,m,k);
                        Qr.fetch(m) += edgeData.phi_yr->get(ell,k)*q.get(i,j,m,k);
                    }
                }

                // Riemann data - aux
                for (int m=1; m<=maux; m++)
                {
                    Auxl.set(m, 0.0 );
                    Auxr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Auxl.fetch(m) += edgeData.phi_yl->get(ell,k)*aux.get(i,j-1,m,k);
                        Auxr.fetch(m) += edgeData.phi_yr->get(ell,k)*aux.get(i,j,m,k);
                    }
                }



                // Solve Riemann problem
                //xedge.set(1, xlower + (double(i)-0.5)*dx );
                //xedge.set(2, ylower + (double(j)-1.0)*dy );
                double s = edgeData.xpts1d->get(ell);
                xedge.set(1, x1p + 0.5*(s+1.0)*(x2p-x1p) );
                xedge.set(2, y1p + 0.5*(s+1.0)*(y2p-y1p) );
               // cout<<"xdiff "<<x1p + 0.5*(1.0)*(x2p-x1p)<<" "<<xlower + (double(i)-0.5)*dx<<endl;
               // cout<<"ydiff "<<y1p + 0.5*(1.0)*(y2p-y1p) <<" "<<ylower + (double(j)-1.0)*dy<<endl;
                

                // TODO: why is this different than what was used earlier in
                // the code? (-DS)
                const double smax_edge = riemannSolver.solve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);

                // TODO: this looks like a potential problem for the upper
                // pragma statement.  See above comment. (-DS)
                //
                smax.set(i, j-1, 2, Max(0.5/mag*smax_edge,smax.get(i, j-1,2) ) );
                smax.set(i, j,   2, Max(0.5/mag*smax_edge,smax.get(i, j,  2) ) );

           //     smax.set(i, j-1,2, Max(abs(smax_edge),smax.get(i, j,2) ));
           //     smax.set(i, j,2,  Max(abs(smax_edge),smax.get(i, j+1,2) ) );

                // Construct fluxes
                for (int m=1; m<=meqn; m++)
                {
                    Gm.set(i,j,  m,ell, Fr.get(m)*mag );
                    Gp.set(i,j-1,m,ell, Fl.get(m)*mag );
                }
            }
          }
    }

    // Compute ``flux differences'' dF and dG    
    const int mlength = Lstar.getsize(3);

#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)    
        for (int j=(2-mbc); j<=(my+mbc-1); j++)
            for (int m=1; m<=mlength; m++)
                for (int k=1; k<=kmax; k++)
                {

               double x1 = xlow + (double(i+1)-1.0)*dx;
               double y1 = ylow + (double(j+1)-1.0)*dy;
               double x2 = xlow + (double(i+1))*dx;
               double y2 = ylow + (double(j+1-1))*dy; 
               double x3 = xlow + (double(i+1))*dx;
               double y3 = ylow + (double(j+1))*dy;
               double x4 = xlow + (double(i+1)-1)*dx;
               double y4 = ylow + (double(j+1))*dy;
               mapc2p(x1,y1);mapc2p(x2,y2);mapc2p(x3,y3);mapc2p(x4,y4);
               double areat=abs(0.5*(x2-x1)*0.5*(y4-y1)-0.5*(x4-x1)*0.5*(y2-y1))*4.0;



                    // 1-direction: dF
                    double F1 = 0.0;
                    double F2 = 0.0;
                    for (int ell=1; ell<=space_order; ell++)
                    {
                        F1 = F1 + edgeData.wght_phi_xr->get(ell,k)*Fm.get(i,j,m,ell);
                        F2 = F2 + edgeData.wght_phi_xl->get(ell,k)*Fp.get(i,j,m,ell);
                    }

                    // 2-direction: dG
                    double G1 = 0.0;
                    double G2 = 0.0;
                    for (int ell=1; ell<=space_order; ell++)
                    {
                        G1 = G1 + edgeData.wght_phi_yr->get(ell,k)*Gm.get(i,j,m,ell);
                        G2 = G2 + edgeData.wght_phi_yl->get(ell,k)*Gp.get(i,j,m,ell);
                    }
                    Lstar.fetch(i,j,m,k) -= (areat)/(dx*dy)*((F2-F1)+ (G2-G1));
                }

    */

Integrator BasisInterator;

#pragma omp parallel
{
        dTensor1 Ql(meqn),Qr(meqn);
        dTensor1 Auxl(maux),Auxr(maux);
        dTensor1 xedge(2);

        RiemannSolver riemannSolver(meqn,maux);
        dTensor1 nvec(2);
        dTensor1 Fl(meqn),Fr(meqn);

      #pragma omp for
    for (int i=(1-mbc); i<=(mx+mbc-1); i++)
    {    
    
        int i1=i+mbc;

        for (int j=(1-mbc); j<=(my+mbc); j++)
        {
        dTensor2 Fr_tmp(meqn,dogParams.get_space_order());
        dTensor2 Fl_tmp(meqn,dogParams.get_space_order());
        int j1=j+mbc;
 
        double xil1 = xlower + (double(i)-1.0)*dx;
        double yil1 = ylower + (double(j)-1.0)*dy;
        double xil2 = xlower + (double(i))*dx;
        double yil2 = ylower + (double(j-1))*dy; 
        double xil3 = xlower + (double(i))*dx;
        double yil3 = ylower + (double(j))*dy;
        double xil4 = xlower + (double(i)-1)*dx;
        double yil4 = ylower + (double(j))*dy;
        vector<double> x=returnright(xil1,yil1,xil2,yil2,xil3,yil3,xil4,yil4,2);
        double x1l = x.at(2);
        double y1l = x.at(3); 
        double x2l = x.at(4);
        double y2l = x.at(5);
        double x3l = x.at(6);
        double y3l = x.at(7);


        double areal=x.at(1);


        double areal2=x.at(0);

        double xir1 = xlower + (double(i+1)-1.0)*dx ;
        double yir1 = ylower + (double(j)-1.0)*dy;
        double xir2 = xlower + (double(i+1))*dx;
        double yir2 = ylower + (double(j-1))*dy; 
        double xir3 = xlower + (double(i+1))*dx;
        double yir3 = ylower + (double(j))*dy;
        double xir4 = xlower + (double(i+1)-1)*dx;
        double yir4 = ylower + (double(j))*dy;
        x=returnleft(xir1,yir1,xir2,yir2,xir3,yir3,xir4,yir4,2);
        double x1r = x.at(2);
        double y1r = x.at(3);
        double x2r = x.at(4);
        double y2r = x.at(5); 
        double x3r = x.at(6);
        double y3r = x.at(7);

        double arear=x.at(0);


        double arear2=x.at(1);
        nvec.set(1, (y1l-y2l) );
        nvec.set(2, (x2l-x1l) ); 

        double length=sqrt(nvec.get(1)*nvec.get(1)+nvec.get(2)*nvec.get(2));
        nvec.set(1, (y1l-y2l)/length );
        nvec.set(2, (x2l-x1l)/length ); 

            for (int ell=1; ell<=space_order; ell++)
            {
                // Riemann data - q
                for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, 0.0 );
                    Qr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {   
                        Ql.fetch(m) += edgeData.phi_xl->get(ell,k)*q.get(i,j,m,k);
                        Qr.fetch(m) += edgeData.phi_xr->get(ell,k)*q.get(i+1,j,m,k);
                    }
                }
                
                // Riemann data - aux
                for (int m=1; m<=maux; m++)
                {
                    Auxl.set(m, 0.0 );
                    Auxr.set(m, 0.0 );
                    Auxl.fetch(m) = aux.get(i,j,m,1);
                    Auxr.fetch(m) = aux.get(i+1,j,m,1);
                   
                }

              double s = edgeData.xpts1d->get(ell);
                xedge.set(1, x2l + 0.5*(s+1.0)*(x1l-x2l) );
                xedge.set(2, y2l + 0.5*(s+1.0)*(y1l-y2l) );
                const double smax_edge
                    = riemannSolver.solve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);
            // Construct fluxes
            for (int m=1; m<=meqn; m++)
            {
                Fr_tmp.set(m,ell, Fr.get(m) );
                Fl_tmp.set(m,ell, Fl.get(m) );
            }

                // TODO: See comments in all ConstructL methods.  This *should* be a 
                // bug when OpenMP loops are turned on. (-DS)
                //
                // If we reverse the role of i and j in the above for loop,
                // then it will be OK, but the tensors don't run as fast when
                // that's the case.
                //
                #pragma omp critical
                {smax.set(i, j,1,  Max(abs(smax_edge)*length, smax.get(i, j,1) ) );}

                #pragma omp critical
                {smax.set(i+1,   j,1,  Max(abs(smax_edge)*length, smax.get(i+1,   j,1) ) );}
                
            }

        for (int m=1; m<=meqn; m++)
            for (int k=1; k<=kmax; k++)
            {
                double Fl_sum = 0.0;
                double Fr_sum = 0.0;
                for (int ell=1; ell<=dogParams.get_space_order(); ell++)
                {   //dont need factor of 2 because of lhs jacobian
                    Fl_sum = Fl_sum + 0.5*edgeData.wgts1d->get(ell)
                        *edgeData.phi_xl->get(ell,k) *Fl_tmp.get(m,ell);
                    Fr_sum = Fr_sum + 0.5*edgeData.wgts1d->get(ell)
                        *edgeData.phi_xr->get(ell,k)*Fr_tmp.get(m,ell);
                }
             //   cout<<" here "<<Fl_sum<<" "<<Fr_sum<<endl;
              //  cout<<"Lstar1= "<<Lstar.get(i,j,m,k)<<" "<<Fl_sum<<endl;
                #pragma critical
                {
                Lstar.set(i,j,m,k,Lstar.get(i,j,m,k)-Fl_sum*length/(abs(areal)+abs(areal2)));           //     cout<<"Lstar2= "<<Lstar.get(i,j,m,k)<<endl;
                Lstar.set(i+1,j,m,k,Lstar.get(i+1,j,m,k)+Fr_sum*length/(abs(arear)+abs(arear2)));
                }
            }

        }
    }
}


    // 2-direction: loop over interior edges and solve Riemann problems

#pragma omp parallel
{
                RiemannSolver riemannSolver(meqn,maux);
        dTensor1 Ql(meqn),Qr(meqn);
        dTensor1 Auxl(maux),Auxr(maux);
        dTensor1 xedge(2);
        dTensor1 nvec(2);
        dTensor1 Fl(meqn),Fr(meqn);
#pragma omp for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
        int i1=i+mbc;
        for (int j=(1-mbc); j<=(my+mbc-1); j++)
        {
        int j1=j+mbc;
        double xil1 = xlower + (double(i)-1.0)*dx;
        double yil1 = ylower + (double(j)-1.0)*dy;
        double xil2 = xlower + (double(i))*dx;
        double yil2 = ylower + (double(j-1))*dy; 
        double xil3 = xlower + (double(i))*dx;
        double yil3 = ylower + (double(j))*dy;
        double xil4 = xlower + (double(i)-1)*dx;
        double yil4 = ylower + (double(j))*dy;

        vector<double> x=returnright(xil1,yil1,xil2,yil2,xil3,yil3,xil4,yil4,1);
        double x1l = x.at(2);
        double y1l = x.at(3); 
        double x2l = x.at(4);
        double y2l = x.at(5);
        double x3l = x.at(6);
        double y3l = x.at(7);
        double areal=x.at(1);
        double areal2=x.at(0);

        double xir1 = xlower + (double(i)-1.0)*dx ;
        double yir1 = ylower + (double(j+1)-1.0)*dy;
        double xir2 = xlower + (double(i))*dx;
        double yir2 = ylower + (double(j+1-1))*dy; 
        double xir3 = xlower + (double(i))*dx;
        double yir3 = ylower + (double(j+1))*dy;
        double xir4 = xlower + (double(i)-1)*dx;
        double yir4 = ylower + (double(j+1))*dy;

        x=returnleft(xir1,yir1,xir2,yir2,xir3,yir3,xir4,yir4,1);
        double x1r = x.at(2);
        double y1r = x.at(3);
        double x2r = x.at(4);
        double y2r = x.at(5); 
        double x3r = x.at(6);
        double y3r = x.at(7);

        double arear=x.at(0);



        double arear2=x.at(1);
        dTensor2 Fr_tmp(meqn,dogParams.get_space_order());
        dTensor2 Fl_tmp(meqn,dogParams.get_space_order());
        nvec.set(1, (y3l-y1l) );
        nvec.set(2, (x1l-x3l) );
        double length=sqrt(nvec.get(1)*nvec.get(1)+nvec.get(2)*nvec.get(2));
        nvec.set(1, (y3l-y1l)/length );
        nvec.set(2, (x1l-x3l)/length ); 



        double G1 = 0.0;
        double G2 = 0.0;
           for (int ell=1; ell<=space_order; ell++)
            {
                // Riemann data - q
                for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, 0.0 );
                    Qr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Ql.fetch(m) += edgeData.phi_yl->get(ell,k)*q.get(i,j,m,k);
                        Qr.fetch(m) += edgeData.phi_yr->get(ell,k)*q.get(i,j+1,m,k);
                    }
                }

                // Riemann data - aux
                for (int m=1; m<=maux; m++)
                {
                    Auxl.set(m, 0.0 );
                    Auxr.set(m, 0.0 );
                    Auxl.fetch(m) = aux.get(i,j,m,1);
                    Auxr.fetch(m) = aux.get(i,j+1,m,1);
                }
               

                // Solve Riemann problem
              double s = edgeData.xpts1d->get(ell);
                xedge.set(1, x1l + 0.5*(s+1.0)*(x3l-x1l) );
                xedge.set(2, y1l + 0.5*(s+1.0)*(y3l-y1l) );

                // TODO: why is this different than what was used earlier in
                // the code? (-DS)
                const double smax_edge = riemannSolver.solve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);

                // TODO: this looks like a potential problem for the upper
                // pragma statement.  See above comment. (-DS)
                //
                #pragma omp critical
                {smax.set(i, j,2, Max(abs(smax_edge)*length,smax.get(i, j,2) ) );}
                #pragma omp critical
                {smax.set(i, j+1,2,  Max(abs(smax_edge)*length,smax.get(i, j+1,2) ) );}
                // Construct fluxes
            for (int m=1; m<=meqn; m++)
            {
                Fr_tmp.set(m,ell, Fr.get(m) );
                Fl_tmp.set(m,ell, Fl.get(m) );
            }

            }

        
        for (int m=1; m<=meqn; m++)
            for (int k=1; k<=kmax; k++)
            {                
                double Gl_sum = 0.0;
                double Gr_sum = 0.0;
                for (int ell=1; ell<=dogParams.get_space_order(); ell++)
                {   //dont need factor of 2 because of lhs jacobian
                    Gl_sum = Gl_sum + 0.5*edgeData.wgts1d->get(ell)
                        *edgeData.phi_yl->get(ell,k) *Fl_tmp.get(m,ell);
                    Gr_sum = Gr_sum + 0.5*edgeData.wgts1d->get(ell)
                        *edgeData.phi_yr->get(ell,k)*Fr_tmp.get(m,ell);
                }
                //cout<<"Lstar4= "<<Lstar.get(i,j,m,k)<<" "<<Gl_sum<<endl;
                #pragma critical
                {
                Lstar.set(i,j,m,k,Lstar.get(i,j,m,k)-Gl_sum*length/(abs(areal)+abs(areal2)));
                Lstar.set(i,j+1,m,k,Lstar.get(i,j+1,m,k)+Gr_sum*length/(abs(arear)+abs(arear2)));   
                }
            }

           }
     }

}

    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part III: compute intra-element contributions
    // ---------------------------------------------------------
    // No need to call this if first-order in space
    if(dogParams.get_space_order()>1)
    {
        L2ProjectGradAdd(1-mbc,mx+mbc,1-mbc,my+mbc,
                space_order,space_order,space_order,space_order,
                &q,&aux,&Lstar,&FluxFunc);
    }
    // ---------------------------------------------------------  


    // ---------------------------------------------------------
    // Part IV: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra(&q,&aux,&Lstar);
    // ---------------------------------------------------------

    if(dogParams.get_space_order()>1)
    {
    BasisInterator.IntegrateBasis(1-mbc,mx+mbc,1-mbc,my+mbc,
                space_order,space_order,space_order,space_order,
                &q,&aux,&Lstar);
    }

    // ---------------------------------------------------------
    // Part V: artificial viscosity limiter
    // ---------------------------------------------------------  
    if (dogParams.get_space_order()>1  &&
            dogParams.using_viscosity_limiter())
    {  ArtificialViscosity(&aux,&q,&Lstar);  }
    // ---------------------------------------------------------
}
