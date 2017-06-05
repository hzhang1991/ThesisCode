#ifndef __DogSolverCart2_h__
#define __DogSolverCart2_h__
#include "DogSolver.h"
#include "DogStateCart2.h"
#include "RKinfo.h"
#include "SDCinfo.h"

class dTensor1;
class dTensor2;
class dTensor3;
class dTensorBC3;
class dTensorBC4;
class dTensorBC5;
class RKinfo;
class DogStateCart2;

class DogSolverCart2: public DogSolverTB
{

    public: // constructor and destructor
        DogSolverCart2(): smax(NULL),L(NULL){};
        virtual ~DogSolverCart2();

    public: // methods
        virtual void init();
        virtual void initParams();

        // these should be consolidated into DogSolve
        //
        // Where's an appropriate spot for a user to implement a new
        // DogSolveUser? (-DS)
        //
        //virtual void DogSolveRK   (double tstart, double tend);
        virtual void DogSolveLxW  (double tstart, double tend);
        virtual void DogSolveSDC  (double tstart, double tend);
        virtual void DogSolveUser (double tstart, double tend);

        // expected to be called between time frames
        virtual void write_restart(int n) const;
        virtual void Output(int n) const;
        virtual void Restart(int nstart);

    public:

        // accessors (fetch means write-access, get means read-access)
        //
        static const DogSolverCart2& get_solver()
        {return (DogSolverCart2&)DogSolver::get_solver();}

        static DogSolverCart2& fetch_solver()
        {return (DogSolverCart2&)DogSolver::fetch_solver();}

        const DogStateCart2& get_state()const
        {return (const DogStateCart2&)DogSolver::get_state();}

        // used by semi_lagrangian DogSolveUser
        DogStateCart2& fetch_state()
        {return (DogStateCart2&)DogSolver::fetch_state();}
        
        dTensorBC3& fetch_smax() {return *smax;}
        const dTensorBC3& get_smax()const{return *smax;}

        dTensorBC4& fetch_L() {return *L;}
        const dTensorBC4& get_L()const{return *L;}

    public:
        // additional functionality for time stepping methods
        RKinfo rk;
        SDCinfo sdc;


    private: // data
        // data that changes but is not necessary to define the state
        // (allocated here for convenience so DogSolve* does not have
        // to reallocate it with each call)
        //
        dTensorBC3* smax;
        dTensorBC4* L;

    protected:
        virtual void initOutputDirectory(int idx, const char* framedir);

    private: // methods
        virtual void initOutputDirectories();
        void advanceTimeStepSDC(double dt,
                dTensorBC4** L,
                dTensorBC4** Lnew,
                dTensorBC4** q,
                dTensorBC5& IL,
                dTensor1& dtvec,
                dTensor1& tvec);
        virtual void saveState();
        virtual void revertToSavedState();
        int advanceTimeStepRK(double dt);

    private: // disabled
        DogSolverCart2& operator=(const DogSolverCart2& in);
        DogSolverCart2(const DogSolverCart2& in);

    private:
        const DogStateCart2& get_state_old()const
        { return (const DogStateCart2&)DogSolver::get_state_old();}
        DogStateCart2& fetch_state_old()
        {return (DogStateCart2&)DogSolver::fetch_state_old();}

    private:
        //virtual void ConstructL(const DogState& state);
        // protected: // uncomment this line if needed
        virtual double GetCFL(double dt) const;
};

void Output(const DogSolverCart2& solver, int n);
void Restart(DogSolverCart2& solver, int nstart);

// public methods
// (could be made static methods of the above class)
//
// core methods
//
void ProjectRightEig(int,const dTensor1&,const dTensor1&,
        const dTensor2&,dTensor2&);
void ProjectLeftEig(int,const dTensor1&,const dTensor1&,
        const dTensor2&,dTensor2&);

// User is responsible for setting boundary data before calling ConstructL:
void ConstructL(
    const dTensorBC4& aux,
    const dTensorBC4& q,
    dTensorBC4& Lstar,
    dTensorBC3& smax);

// ----------------------
// user-defined callbacks
// ----------------------
//
// limiter
//
void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
        void (*ProjectRightEig)(int,const dTensor1&,
            const dTensor1&,const dTensor2&, dTensor2&),
        void (*ProjectLeftEig)(int,const dTensor1&,
            const dTensor1&,const dTensor2&, dTensor2&));

//
// event hooks
//
void AfterQinit(DogSolverCart2& solver);

///////////////////////////////////////////////////////////////////////////////
// EVENTUALLY only the DogSolverCart2 argument should remain and the rest
// should be deleted. In order to do this, however, need to make some
// changes to DogSolveRK and DogSolveSDC
///////////////////////////////////////////////////////////////////////////////
void BeforeStep  (double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver);
void AfterStep   (double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver);
void AfterReject (double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver);
//////
void AfterUpdateSoln (const dTensorBC4& aux, dTensorBC4& q,double dt, double beta);
void AfterFullTimeStep(DogSolverCart2& solver);
void ConSoln         (const dTensorBC4& aux, const dTensorBC4& q, double t);

//
// called on startup
//
void InitApp();
double Qinit_restart(int nstart, dTensorBC4& q);

#endif
