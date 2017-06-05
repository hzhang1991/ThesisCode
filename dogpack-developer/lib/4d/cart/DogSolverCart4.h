#ifndef __DogSolverCart4_h__
#define __DogSolverCart4_h__
#include "DogSolver.h"
#include "DogStateCart4.h"
#include "RKinfo.h"

class dTensor1;
class dTensor2;
class dTensor3;
class dTensorBC5;   // for smax
class dTensorBC6;   // State variables
class RKinfo;
class DogStateCart4;

class DogSolverCart4: public DogSolverTB
{

    public: // constructor and destructor
        DogSolverCart4(): smax(NULL),L(NULL){};
        virtual ~DogSolverCart4();

    public: // methods
        virtual void init();
        virtual void initParams();

        // these should be consolidated into DogSolve
        //
        // Where's an appropriate spot for a user to implement a new
        // DogSolveUser? (-DS)
        //
        //virtual void DogSolveRK   (double tstart, double tend);
        virtual void DogSolveSDC  (double tstart, double tend);
        virtual void DogSolveLxW  (double tstart, double tend);
        virtual void DogSolveUser (double tstart, double tend);

        // expected to be called between time frames
        virtual void write_restart(int n) const;
        virtual void Output(int n) const;
        virtual void Restart(int nstart);

    public:

        // accessors (fetch means write-access, get means read-access)
        //
        static const DogSolverCart4& get_solver()
        {return (DogSolverCart4&)DogSolver::get_solver();}

        static DogSolverCart4& fetch_solver()
        {return (DogSolverCart4&)DogSolver::fetch_solver();}

        const DogStateCart4& get_state()const
        {return (const DogStateCart4&)DogSolver::get_state();}

        // used by semi_lagrangian DogSolveUser
        DogStateCart4& fetch_state()
        {return (DogStateCart4&)DogSolver::fetch_state();}

        dTensorBC5& fetch_smax() {return *smax;}
        const dTensorBC5& get_smax()const{return *smax;}

        dTensorBC6& fetch_L() {return *L;}
        const dTensorBC6& get_L()const{return *L;}

    public:
        // additional functionality for time stepping methods
        RKinfo rk;        


    private: // data
        // data that changes but is not necessary to define the state
        // (allocated here for convenience so DogSolve* does not have
        // to reallocate it with each call)
        //
        dTensorBC5* smax;
        dTensorBC6* L;

    protected:
        virtual void initOutputDirectory(int idx, const char* framedir);

    private: // methods
        virtual void initOutputDirectories();
        virtual void saveState();
        virtual void revertToSavedState();
        int advanceTimeStepRK(double dt);

    private: // disabled
        DogSolverCart4& operator=(const DogSolverCart4& in);
        DogSolverCart4(const DogSolverCart4& in);

    private:
        const DogStateCart4& get_state_old()const
        { return (const DogStateCart4&)DogSolver::get_state_old();}
        DogStateCart4& fetch_state_old()
        {return (DogStateCart4&)DogSolver::fetch_state_old();}

    private:
        //virtual void ConstructL(const DogState& state);
        // protected: // uncomment this line if needed
        virtual double GetCFL(double dt) const;
};

void Output(const DogSolverCart4& solver, int n);
void Restart(DogSolverCart4& solver, int nstart);

// public methods
// (could be made static methods of the above class)
//
// core methods
//
void ProjectRightEig(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&);
void ProjectLeftEig(int,
        const dTensor1&,
        const dTensor1&,
        const dTensor2&,
        dTensor2&);

// User is responsible for setting boundary data before calling ConstructL:
void ConstructL(const dTensorBC6& aux, const dTensorBC6& q, dTensorBC6& Lstar,
        dTensorBC5& smax);

// ----------------------
// user-defined callbacks
// ----------------------
//
// limiter
//
void ApplyLimiter(dTensorBC6& aux, dTensorBC6& q,
        void (*ProjectRightEig)(int,
            const dTensor1&,
            const dTensor1&,
            const dTensor2&, 
            dTensor2&),
        void (*ProjectLeftEig)(int,
            const dTensor1&,
            const dTensor1&,
            const dTensor2&, 
            dTensor2&));

//
// event hooks
//
void AfterQinit(DogSolverCart4& solver);

///////////////////////////////////////////////////////////////////////////////
// EVENTUALLY only the DogSolverCart4 argument should remain and the rest
// should be deleted. In order to do this, however, need to make some
// changes to DogSolveRK and DogSolveSDC
///////////////////////////////////////////////////////////////////////////////
void BeforeStep         (double dt, DogSolverCart4& solver);
void AfterStep          (double dt, DogSolverCart4& solver);
void AfterReject        (double dt, DogSolverCart4& solver);
void AfterUpdateSoln    (double dt, DogSolverCart4& solver);
void ConSoln            (double t);

//////
void AfterFullTimeStep(DogSolverCart4& solver);

//
// called on startup
//
void InitApp();
double Qinit_restart(int nstart, dTensorBC5& q);

#endif
