#ifndef _LEGENDRE2D_H_
#define _LEGENDRE2D_H_

#include<cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Interval.h"
#include "dog_math.h"
#include "edge_data.h"
#include "Quadrature.h"

class dTensor1;
class dTensor2;
class dTensor3;
class dTensorBC4;
class IntervalArray;
class edge_data;

void SetQuadrature(const int mpoints1d, dTensor1& wgt, dTensor2& spts);
void SetLegendreExtrema(dTensor1& phi);
void SetLegendreAtPoints(const dTensor2& spts, dTensor2& phi);
void SetLegendre(const int mpoints1d,
        const dTensor2& spts,
        dTensor2& phi);
void MP_SetLegendre(const int mpoints1d,
        const dTensor2& spts,
        dTensor2& phi);
void SetLegendreGrad(const double dx, const double dy,
        const dTensor2& spts, dTensor2& phi_x, dTensor2& phi_y);
void MP_SetLegendreGrad(const double dx, const double dy,
        const dTensor2& spts, dTensor2& phi_x, dTensor2& phi_y);
void project_onto_locally_divergence_free_subspace(dTensorBC4& q);

// singleton
//
class Legendre2d
{

    private:
        edge_data* edgeData; // used for integrating along edges
        //dTensor2* gaussLobattoSpts;
        //dTensor1* legendreExtrema;
        dTensor1* quadratureWeights; //wgt
        dTensor2* quadraturePoints; //spts
        dTensor2* phi; // values of legendre basis at quadrature points
        dTensor2* wgt_phi_transpose ; // quadratureWeights*quadraturePoints^T
        dTensor3* phi_divfree; // values of divfree basis at quadrature points
        dTensor1* phi_min_ptr; // minimum of legendre basis functions
        dTensor1* phi_max_ptr; // maximum of legendre basis functions
        IntervalArray* phi_interval; // array of intervals [phi_min(k), phi_max(k)]
        //dTensor1* phi_maxabs; // infinity-norm of legendre basis functions
        //
        // for one-order-lower integration
        // (used when integrating gradient in L2ProjectGrad)
        //
        dTensor1* quadratureWeightsCoarse; //wgt
        dTensor2* quadraturePointsCoarse; //spts
        dTensor2* phiCoarse; // legendre values at quadrature points
        dTensor2* phi_x; // legendre_x_values (x-derivatives)
        dTensor2* phi_y; // legendre_y_values (y-derivatives)
        //
        // points at which to maintain positivity
        // (a sufficiently rich set with a sufficiently short
        // time step guarantees that positivity of the
        // average is maintained.)
        dTensor2* positivityPoints;
        int numPositivityPoints;
        dTensor2* phiAtPositivityPoints; // legendre values at positivity points
        dTensor2* riemannPoints;
        int numRiemannPoints;
        dTensor2* phiAtRiemannPoints; // legendre values at positivity points
    private:
        ~Legendre2d();
        Legendre2d(); // disable constructor of singleton
        void setPositivityPoints();
        void setRiemannPoints();
        void SetLegendreIntervals();
    public:

        void check_edgeData() const;

        static const dTensor1& get_quadratureWeights() //const
        { return *(instance().quadratureWeights);}

        static const dTensor2& get_quadraturePoints() //const
        { return *(instance().quadraturePoints);}

        static const dTensor2& get_wgt_phi_transpose() //const
        { return *(instance().wgt_phi_transpose);}

        static const dTensor2& get_phi() //const
        {return *(instance().phi);}

        static const dTensor3& get_phi_divfree() //const
        {return *(instance().phi_divfree);}
        //
        //static int get_numCoarsePoints() //const
        //  {return instance().numCoarsePoints;}
        static const dTensor2& get_phiCoarse() //const
        {return *(instance().phiCoarse);}

        static int get_numPositivityPoints() //const
        {return instance().numPositivityPoints;}

        static const dTensor2& get_phiAtPositivityPoints() //const
        {return *(instance().phiAtPositivityPoints);}

        static const dTensor2& get_positivityPoints() //const // used?
        {return *(instance().positivityPoints);}

        static const dTensor2& get_phiAtRiemannPoints() //const
        {return *(instance().phiAtRiemannPoints);}

        static const dTensor2& get_riemannPoints() //const // used?
        {return *(instance().riemannPoints);}

        static const edge_data& get_edgeData() // const
        {return *(instance().edgeData);}

        static const IntervalArray& get_phi_interval() //const
        {return *(instance().phi_interval);}

        static Legendre2d& instance(); //const;
};

#endif
