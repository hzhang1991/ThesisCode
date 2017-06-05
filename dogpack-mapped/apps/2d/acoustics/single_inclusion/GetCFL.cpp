#include "tensors.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogSolverCart2.h"
#include <iostream>
using namespace std;



void mapc2p(double& xc,double& yc);
vector<double> returnleft(double x1l,double y1l,double x2l,double y2l,double x3l,double y3l,double x4l,double y4l,int);
vector<double> returnright(double x1l,double y1l,double x2l,double y2l,double x3l,double y3l,double x4l,double y4l,int);


double DogSolverCart2::GetCFL(double dt) const
{
const dTensorBC2& smax = get_smax();

    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();
    const double xlower = dogParamsCart2.get_xlow();
    const double ylower = dogParamsCart2.get_ylow();
    int mx=dogParamsCart2.get_mx();
    int my=dogParamsCart2.get_my();
    const int mbc = dogParamsCart2.get_mbc();
    double cfl=0.0;



#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    {    

        for (int j=(1-mbc); j<=(my+mbc); j++)
        {


        double xil1 = xlower + (double(i)-1.0)*dx;
        double yil1 = ylower + (double(j)-1.0)*dy;
        double xil2 = xlower + (double(i))*dx;
        double yil2 = ylower + (double(j-1))*dy; 
        double xil3 = xlower + (double(i))*dx;
        double yil3 = ylower + (double(j))*dy;
        double xil4 = xlower + (double(i-1))*dx;
        double yil4 = ylower + (double(j))*dy;

        vector<double> x=returnleft(xil1,yil1,xil2,yil2,xil3,yil3,xil4,yil4,1);

        //these components of x vector contain the areas of the two triangles that make up half of the quad
        double areal=x.at(1);


        double areal2=x.at(0);

        double tmp = smax.get(i,j);

       double area=min(abs(x.at(0)),abs(x.at(1)));
       cfl = Max(0.5*dt*tmp/area, cfl);

        }
    }

    return cfl;
}
