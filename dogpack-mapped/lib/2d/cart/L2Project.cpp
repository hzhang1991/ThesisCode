#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "L2ProjectInline2d.h"
#include "Legendre2d.h"
#include <iostream>
#include <fstream>

// -------------------------------------------------------------
// Routine for computing the L2-projection of an input function
// onto an orthonormal Legendre basis
// -------------------------------------------------------------

void mapc2p(double& xc,double& yc);

static void L2ProjectAdd(const bool add,
        const int istart,
        const int iend,
        const int jstart,
        const int jend,
        const int QuadOrder,
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,    
        const dTensorBC4* qin,
        const dTensorBC4* auxin,
        dTensorBC4* fout,
        void (*Func)(),
        bool has_data, 
        void* data)
{
    if(!has_data) assert(!data);

    // dx and dy
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double xlow = dogParamsCart2.get_xlow();
    const double ylow = dogParamsCart2.get_ylow();
    const double xmax = dogParamsCart2.get_xhigh();
    const double ymax = dogParamsCart2.get_yhigh();
    // mbc
    const int mbc = qin->getmbc();
    assert_eq(mbc,auxin->getmbc());
    assert_eq(mbc,fout->getmbc());

    // qin variable
    const int       mx = qin->getsize(1);
    const int       my = qin->getsize(2);
    const int     meqn = qin->getsize(3);
    const int kmax_qin = qin->getsize(4);
    assert_eq(kmax_qin,(BasisOrder_qin*BasisOrder_qin));

    // auxin variable
    assert_eq(mx,auxin->getsize(1));
    assert_eq(my,auxin->getsize(2));
    const int       maux = auxin->getsize(3);
    const int kmax_auxin = auxin->getsize(4);
    assert_eq(kmax_auxin,(BasisOrder_auxin*BasisOrder_auxin));

    // fout variables
    //  TODO - why assume this has the same size?  What if we want to only
    //  project onto a single line, e.g. for padding boundary cell data? (-DS)
    // assert_eq(mx,fout->getsize(1));  
    assert_eq(my,fout->getsize(2));
    const int mcomps_out = fout->getsize(3);
    const int  kmax_fout = fout->getsize(4);
    assert_eq(kmax_fout,(BasisOrder_fout*BasisOrder_fout));

    // starting and ending indeces
    assert_ge(istart,1-mbc);
    assert_le(iend,mx+mbc);
    assert_ge(jstart,1-mbc);
    assert_le(jend,my+mbc);

    // number of quadrature points
    assert_ge(QuadOrder,1);
    assert_le(QuadOrder,20);
    int QuadOrder_MOD = QuadOrder+1;
    switch(QuadOrder)
    {
        case 7:
            QuadOrder_MOD = 8;
            break;
        case 9:
            QuadOrder_MOD = 10;
            break;
        case 11:
            QuadOrder_MOD = 12;
            break;
        case 13:
            QuadOrder_MOD = 14;
            break;
        case 15:
            QuadOrder_MOD = 16;
            break;
        case 17:
            QuadOrder_MOD = 18;
            break;
        case 19:
            QuadOrder_MOD = 20;
            break;
    }
    const int mpoints = QuadOrder_MOD*QuadOrder_MOD;

    // set quadrature point and weight information
    void SetQuadWgtsPts(const int, dTensor1&, dTensor2&);
    dTensor1 wgt(mpoints);
    dTensor2 spts(mpoints, 2);
    SetQuadWgtsPts(QuadOrder_MOD, wgt, spts);

    // Loop over each quadrature point to construct Legendre polys
    const int kmax = iMax(iMax(kmax_qin,kmax_auxin),kmax_fout);
    dTensor2 phi(mpoints, kmax);
    void SetLegendrePolys(const int, const int, const dTensor2&, dTensor2&);
    SetLegendrePolys(mpoints, kmax, spts, phi);

    // For efficiency compute weight*phi and then take transpose
    dTensor2 wgt_phi_transpose(kmax,mpoints);
    for(int mp=1;mp<=mpoints;mp++)
    for(int k=1;k<=kmax;k++)
    {
        wgt_phi_transpose.set(k,mp, wgt.get(mp)*phi.get(mp,k) );
    }

    // ------------------------------------------------------------- //
    // Loop over every grid cell indexed by user supplied parameters //
    // described by istart...iend, jstart...jend                     // 
    // ------------------------------------------------------------- //
ofstream myfile;
myfile.open ("mapped.txt");
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {

        // Local storage for q, aux and xpts (all passed into user supplied
        // function)
        dTensor2   qvals(mpoints, meqn);
        dTensor2 auxvals(mpoints, maux);
        dTensor2    xpts(mpoints, 2);

        // Flux function and its Jacobian:
        dTensor2   fvals(mpoints, mcomps_out);

        for (int j=jstart; j<=jend; j++)
        {
            //find center of current cell
            const double xc = dogParamsCart2.get_xc(i);
            const double yc = dogParamsCart2.get_yc(j);

        const double xce = xlow + (double(i)-0.5)*dx;
        const double yce = ylow + (double(j)-0.5)*dy;

            int jm=j;
            double xi1 = xlow+(i-1)*dx;
            double eta1= ylow+(jm-1)*dy;
            double xi2 = xlow+(i)*dx;
            double eta2= ylow+(jm-1)*dy;
            double xi3 = xlow+(i)*dx;
            double eta3= ylow+(jm)*dy;
            double xi4 = xlow+(i-1)*dx;
            double eta4= ylow+(jm)*dy;

            double xp1 = xlow+(i-1)*dx;
            double yp1= ylow+(jm-1)*dy;
            double xp2 = xlow+(i)*dx;
            double yp2= ylow+(jm-1)*dy;
            double xp3 = xlow+(i)*dx;
            double yp3= ylow+(jm)*dy;
            double xp4 = xlow+(i-1)*dx;
            double yp4= ylow+(jm)*dy;

            mapc2p(xp1,yp1);mapc2p(xp2,yp2);mapc2p(xp3,yp3);mapc2p(xp4,yp4);
            if(xi1>=xlow && eta1>=ylow && xi1<=xmax && eta1<=ymax)
            myfile << xp1<<" "<<yp1<<"\n";
            double xc1=(xp1+xp2+xp3+xp4)/4.0;
            double yc0=(yp1+yp2+yp3+yp4)/4.0;
            double xmax=max(max(xp1,xp2),max(xp3,xp4));double xmin=min(min(xp1,xp2),min(xp3,xp4));
            double ymax=max(max(yp1,yp2),max(yp3,yp4));double ymin=min(min(yp1,yp2),min(yp3,yp4));
            double dx1=xmax-xmin;double dy1=ymax-ymin;


            // Compute q, aux and fvals at each Gaussian quadrature point
            // for this current cell indexed by (i,j)
            // Save results into dTensor2 qvals, auxvals and fvals.
    void mapc2p(double& xc,double& yc);
    double mapx(double xi,double eta,double xp1,double xp2,double xp3,double xp4,double yp1,double yp2,double yp3,double yp4);
    double mapy(double xi,double eta,double xp1,double xp2,double xp3,double xp4,double yp1,double yp2,double yp3,double yp4);
            L2ProjectInline2d::set_vals_at_each_Gaussian_quadrature_point(i, j, 
                    mpoints, meqn, maux, 
                    kmax_qin, kmax_auxin, 
                    xce, yce, dx, dy, 
                    spts, phi, qin, auxin, 
                    xpts, qvals, auxvals);
            // Evaluate Func at Gaussian quadrature point
            if(has_data)
            {
                ((void (*)(const dTensor2&, const dTensor2&, 
                           const dTensor2&, dTensor2&,  void* data))
                 Func)(xpts, qvals, auxvals, fvals, data);
            }
            else
            {
                ((void (*)(const dTensor2&, const dTensor2&, 
                           const dTensor2&, dTensor2&))
                 Func)(xpts, qvals, auxvals, fvals);
            }

            // Evaluate integral on current cell (project onto Legendre basis)
            // using Gaussian Quadrature for the integration
            L2ProjectInline2d::integrate_on_current_cell(add, i, j, 
                    mcomps_out, kmax_fout, mpoints, 
                    wgt_phi_transpose, 
                    fvals, fout);
        }
    }

myfile.close();
}

// L2Project
void L2Project(const int istart,
        const int iend,
        const int jstart,
        const int jend,
        const int QuadOrder,
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,    
        const dTensorBC4* qin,
        const dTensorBC4* auxin,
        dTensorBC4* fout,
        void (*Func)(const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor2&))
{
    L2ProjectAdd(false,istart,iend,jstart,jend,
            QuadOrder,BasisOrder_qin,BasisOrder_auxin,
            BasisOrder_fout,qin,auxin,fout,
            ((void(*)()) Func),false,0);
}
void L2Project(int istart, int iend, int jstart, int jend,
        const dTensorBC4& q,
        const dTensorBC4& aux, dTensorBC4& Fout,
        void (*Func)(const dTensor2& xpts,
            const dTensor2& qvals,
            const dTensor2& auxvals,
            dTensor2& source))
{
    const int space_order = dogParams.get_space_order();
    L2ProjectAdd(false,istart,iend,jstart,jend,
            space_order,space_order,space_order,space_order,
            &q,&aux,&Fout,
            ((void(*)()) Func),false,0);
}

// L2Project version to support a callback function that
// receives extra data registered by the user.
//
// This mechanism could be used to pass DogSolver::solver if
// we someday decide that it should no longer be a singleton;
// for now it is accessible via DogSolver::get_solver();
//
// This mechanism can be used to pass time increment
// information. Note that the time increment is available
// via DogSolver::get_solver().get_dt(), and the time of the
// current stage is set in DogSolver::advanceTimeStageRK() and
// is available via a call to DogSolver::get_time_hack().
//
void L2Project_extra(const int istart,
        const int iend,
        const int jstart,
        const int jend,
        const int QuadOrder,
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,    
        const dTensorBC4* qin,
        const dTensorBC4* auxin,
        dTensorBC4* fout,
        void (*Func)(const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor2&, void* data),
        void* data)
{
    L2ProjectAdd(false,istart,iend,jstart,jend,
            QuadOrder,BasisOrder_qin,BasisOrder_auxin,
            BasisOrder_fout,qin,auxin,fout,
            ((void(*)()) Func),true,data);
}
void L2Project_extra(int istart, int iend, int jstart, int jend,
        const dTensorBC4& q,
        const dTensorBC4& aux, dTensorBC4& Fout,
        void (*Func)(const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor2&, void* data),
        void* data)
{
    const int space_order = dogParams.get_space_order();
    L2ProjectAdd(false,istart,iend,jstart,jend,
            space_order,space_order,space_order,space_order,
            &q,&aux,&Fout,
            ((void(*)()) Func),true,data);
}

// Set quadrature weights and points
void SetQuadWgtsPts(const int mpoints1d, dTensor1& wgt, dTensor2& spts)
{
    {
        const int mpoints = mpoints1d*mpoints1d;
        assert_eq(mpoints, wgt.getsize());
        assert_eq(mpoints, spts.getsize(1));
    }

    const int md2 = mpoints1d/2;
    dTensor1 w1d(mpoints1d),x1d(mpoints1d);

    switch (mpoints1d)
    {
        case 1:
            w1d.set(1, 2.0 );

            x1d.set(1, 0.0 );
            break;

        case 2:
            w1d.set(1,  1.0 );
            w1d.set(2,  1.0 );

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );
            break;

        case 3:
            w1d.set(1,  5.0/9.0 );
            w1d.set(2,  8.0/9.0 );
            w1d.set(3,  5.0/9.0 );

            x1d.set(1,  sq3/sq5 );
            x1d.set(2,  0.0 );
            x1d.set(3, -sq3/sq5 );
            break;

        case 4:
            w1d.set(1, (18.0 - sq3*sq10)/36.0 );
            w1d.set(2, (18.0 + sq3*sq10)/36.0 );
            w1d.set(3, w1d.get(2) );
            w1d.set(4, w1d.get(1) );

            x1d.set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            x1d.set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            x1d.set(3, -x1d.get(2) );
            x1d.set(4, -x1d.get(1) );       
            break;

        case 5:      
            w1d.set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            w1d.set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            w1d.set(3, 128.0/225.0 );
            w1d.set(4, w1d.get(2) );
            w1d.set(5, w1d.get(1) );

            x1d.set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            x1d.set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );
            break;

        case 6:
            w1d.set(1, 0.1713244923791703450402961 );
            w1d.set(2, 0.3607615730481386075698335 );
            w1d.set(3, 0.4679139345726910473898703 );      

            x1d.set(1, 0.9324695142031520278123016 );
            x1d.set(2, 0.6612093864662645136613996 );
            x1d.set(3, 0.2386191860831969086305017 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 8:
            w1d.set(1, 0.1012285362903762591525314 );
            w1d.set(2, 0.2223810344533744705443560 );
            w1d.set(3, 0.3137066458778872873379622 );
            w1d.set(4, 0.3626837833783619829651504 );

            x1d.set(1, 0.9602898564975362316835609 );
            x1d.set(2, 0.7966664774136267395915539 );
            x1d.set(3, 0.5255324099163289858177390 );
            x1d.set(4, 0.1834346424956498049394761 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 10:
            w1d.set(1, 0.0666713443086881375935688 );
            w1d.set(2, 0.1494513491505805931457763 );
            w1d.set(3, 0.2190863625159820439955349 );
            w1d.set(4, 0.2692667193099963550912269 );
            w1d.set(5, 0.2955242247147528701738930 );

            x1d.set(1, 0.9739065285171717200779640 );
            x1d.set(2, 0.8650633666889845107320967 );
            x1d.set(3, 0.6794095682990244062343274 );
            x1d.set(4, 0.4333953941292471907992659 );
            x1d.set(5, 0.1488743389816312108848260 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 12:
            w1d.set(1, 0.0471753363865118271946160 );
            w1d.set(2, 0.1069393259953184309602547 );
            w1d.set(3, 0.1600783285433462263346525 );
            w1d.set(4, 0.2031674267230659217490645 );
            w1d.set(5, 0.2334925365383548087608499 );
            w1d.set(6, 0.2491470458134027850005624 );

            x1d.set(1, 0.9815606342467192506905491 );
            x1d.set(2, 0.9041172563704748566784659 );
            x1d.set(3, 0.7699026741943046870368938 );
            x1d.set(4, 0.5873179542866174472967024 );
            x1d.set(5, 0.3678314989981801937526915 );
            x1d.set(6, 0.1252334085114689154724414 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 14:
            w1d.set(1, 0.0351194603317518630318329 );
            w1d.set(2, 0.0801580871597602098056333 );
            w1d.set(3, 0.1215185706879031846894148 );
            w1d.set(4, 0.1572031671581935345696019 );
            w1d.set(5, 0.1855383974779378137417166 );
            w1d.set(6, 0.2051984637212956039659241 );
            w1d.set(7, 0.2152638534631577901958764 );

            x1d.set(1, 0.9862838086968123388415973 );
            x1d.set(2, 0.9284348836635735173363911 );
            x1d.set(3, 0.8272013150697649931897947 );
            x1d.set(4, 0.6872929048116854701480198 );
            x1d.set(5, 0.5152486363581540919652907 );
            x1d.set(6, 0.3191123689278897604356718 );
            x1d.set(7, 0.1080549487073436620662447 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 16:
            w1d.set(1, 0.0271524594117540948517806 );
            w1d.set(2, 0.0622535239386478928628438 );
            w1d.set(3, 0.0951585116824927848099251 );
            w1d.set(4, 0.1246289712555338720524763 );
            w1d.set(5, 0.1495959888165767320815017 );
            w1d.set(6, 0.1691565193950025381893121 );
            w1d.set(7, 0.1826034150449235888667637 );
            w1d.set(8, 0.1894506104550684962853967 );

            x1d.set(1, 0.9894009349916499325961542 );
            x1d.set(2, 0.9445750230732325760779884 );
            x1d.set(3, 0.8656312023878317438804679 );
            x1d.set(4, 0.7554044083550030338951012 );
            x1d.set(5, 0.6178762444026437484466718 );
            x1d.set(6, 0.4580167776572273863424194 );
            x1d.set(7, 0.2816035507792589132304605 );
            x1d.set(8, 0.0950125098376374401853193 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 18:
            w1d.set(1, 0.0216160135264833103133427 );
            w1d.set(2, 0.0497145488949697964533349 );
            w1d.set(3, 0.0764257302548890565291297 );
            w1d.set(4, 0.1009420441062871655628140 );
            w1d.set(5, 0.1225552067114784601845191 );
            w1d.set(6, 0.1406429146706506512047313 );
            w1d.set(7, 0.1546846751262652449254180 );
            w1d.set(8, 0.1642764837458327229860538 );
            w1d.set(9, 0.1691423829631435918406565 );

            x1d.set(1, 0.9915651684209309467300160 );
            x1d.set(2, 0.9558239495713977551811959 );
            x1d.set(3, 0.8926024664975557392060606 );
            x1d.set(4, 0.8037049589725231156824175 );
            x1d.set(5, 0.6916870430603532078748911 );
            x1d.set(6, 0.5597708310739475346078715 );
            x1d.set(7, 0.4117511614628426460359318 );
            x1d.set(8, 0.2518862256915055095889729 );
            x1d.set(9, 0.0847750130417353012422619 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        case 20:
            w1d.set(1,  0.0176140071391521183118620 );
            w1d.set(2,  0.0406014298003869413310400 );
            w1d.set(3,  0.0626720483341090635695065 );
            w1d.set(4,  0.0832767415767047487247581 );
            w1d.set(5,  0.1019301198172404350367501 );
            w1d.set(6,  0.1181945319615184173123774 );
            w1d.set(7,  0.1316886384491766268984945 );
            w1d.set(8,  0.1420961093183820513292983 );
            w1d.set(9,  0.1491729864726037467878287 );
            w1d.set(10, 0.1527533871307258506980843 );

            x1d.set(1,  0.9931285991850949247861224 );
            x1d.set(2,  0.9639719272779137912676661 );
            x1d.set(3,  0.9122344282513259058677524 );
            x1d.set(4,  0.8391169718222188233945291 );
            x1d.set(5,  0.7463319064601507926143051 );
            x1d.set(6,  0.6360536807265150254528367 );
            x1d.set(7,  0.5108670019508270980043641 );
            x1d.set(8,  0.3737060887154195606725482 );
            x1d.set(9,  0.2277858511416450780804962 );
            x1d.set(10, 0.0765265211334973337546404 );

            for (int k=1; k<=md2; k++)
            { w1d.set(md2+k,  w1d.get(md2+1-k) ); }      
            for (int k=1; k<=md2; k++)
            { x1d.set(md2+k, -x1d.get(md2+1-k) ); }
            break;

        default:
            unsupported_value_error(mpoints1d);
    }

    // Tensor product Gaussian Quadrature
    int k=0;
    for (int m1=1; m1<=(mpoints1d); m1++)
        for (int m2=1; m2<=(mpoints1d); m2++)
        {
            k = k+1;
            wgt.set(k,  w1d.get(m1)*w1d.get(m2) );

            spts.set(k,1, x1d.get(m1) );
            spts.set(k,2, x1d.get(m2) );
        }
}


// Loop over each quadrature point to construct Legendre polys
void SetLegendrePolys(const int mpoints, const int kmax, const dTensor2& spts, dTensor2& phi)
{
    assert_eq(mpoints,spts.getsize(1));
    assert_eq(mpoints,phi.getsize(1));
    assert_eq(kmax,phi.getsize(2));

    for (int m=1; m<=mpoints; m++)
    {
        // grid point (x,y)
        const double xi   = spts.get(m,1);
        const double eta  = spts.get(m,2);
        const double xi2  = xi*xi;
        const double xi3  = xi*xi2;
        const double xi4  = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;             

        // Legendre basis functions at each gaussian quadrature point in the
        // interval [-1,1]x[-1,1].
        switch( kmax )
        {
            case 25:  // fifth order                                 
                phi.set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
                phi.set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq7*(2.5*xi3 - 1.5*xi) );
                phi.set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq5*(1.5*xi2 - 0.5) );
                phi.set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq5*(1.5*eta2 - 0.5) );
                phi.set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*xi );
                phi.set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*eta );
                phi.set( m,18, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi.set( m,17, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );

            case 16:  // fourth order
                phi.set( m,16, sq7*sq7*(2.5*eta3 - 1.5*eta)*(2.5*xi3 - 1.5*xi) );
                phi.set( m,15, sq5*sq7*(2.5*eta3 - 1.5*eta)*(1.5*xi2 - 0.5) );
                phi.set( m,14, sq5*sq7*(2.5*xi3 - 1.5*xi)*(1.5*eta2 - 0.5) );
                phi.set( m,13, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                phi.set( m,12, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );
                phi.set( m,11, sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( m,10,  sq7*(2.5*xi3 - 1.5*xi) );

            case 9:  // third order
                phi.set( m,9, 5.0*(1.5*xi2 - 0.5)*(1.5*eta2 - 0.5) );
                phi.set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi.set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );
                phi.set( m,6,  sq5*(1.5*eta2 - 0.5) );
                phi.set( m,5,  sq5*(1.5*xi2 - 0.5) );

            case 4:  // second order                
                phi.set( m,4,  3.0*xi*eta );                  
                phi.set( m,3, sq3*eta );
                phi.set( m,2, sq3*xi  );

            case 1:  // first order
                phi.set( m,1, 1.0 );

                break;

            default:
                unsupported_value_error(kmax);
        }

    }

}
