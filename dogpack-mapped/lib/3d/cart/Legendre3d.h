#ifndef _LEGENDRE3D_H_
#define _LEGENDRE3D_H_

#include<cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart3.h"
#include "dog_math.h"
#include "FaceData.h"
#include "Quadrature.h"   // 1D quadrature rules

class dTensor1;
class dTensor2;
class dTensor3;
class dTensorBC5;
class FaceData;

// -----------------------------------------------------------------------
// FUNCTIONS USED TO SET QUADRATURE POINTS AND WEIGHTS, AS WELL AS 
// LEGENDRE FUNCTIONS AND DERIVATIVES AT QUADRATURE POINTS
// -----------------------------------------------------------------------
void SetQuadWgtsPts(const int mpoints1d, 
		    dTensor1& wgt, 
		    dTensor2& spts);
void SetLegendrePolys(const int mpoints, 
		      const int kmax, 
		      const dTensor2& spts,
		      dTensor2& phi);
void SetLegendrePolysGrad(const double dx,
			  const double dy,
			  const double dz,
			  const int mpoints, 
			  const int kmax, 
			  const dTensor2& spts, 
			  dTensor2& phi_x,
			  dTensor2& phi_y,
			  dTensor2& phi_z);

// -----------------------
// Legendre3d singleton
// -----------------------
class Legendre3d
{

    private:
        FaceData* faceData; // used for integrating along edges
        ~Legendre3d();
        Legendre3d(); // disable constructor of singleton

    public:
        static const FaceData& get_faceData() // const
        {return *(instance().faceData);}

        static Legendre3d& instance(); //const;
};

#endif
