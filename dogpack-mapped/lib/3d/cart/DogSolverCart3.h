#ifndef __DogSolverCart3_h__
#define __DogSolverCart3_h__
#include "DogSolver.h"
#include "DogStateCart3.h"
#include "RKinfo.h"

class dTensor1;
class dTensor2;
class dTensor3;
class dTensorBC4;
class dTensorBC5;
class RKinfo;
class DogStateCart3;

class DogSolverCart3: public DogSolverTB
{

    public: // constructor and destructor
        DogSolverCart3(): smax(NULL),L(NULL){};
        virtual ~DogSolverCart3();

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
        static const DogSolverCart3& get_solver()
        {return (DogSolverCart3&)DogSolver::get_solver();}

        static DogSolverCart3& fetch_solver()
        {return (DogSolverCart3&)DogSolver::fetch_solver();}

        const DogStateCart3& get_state()const
        {return (const DogStateCart3&)DogSolver::get_state();}

        // used by semi_lagrangian DogSolveUser
        DogStateCart3& fetch_state()
        {return (DogStateCart3&)DogSolver::fetch_state();}
        
        dTensorBC4& fetch_smax() {return *smax;}
        const dTensorBC4& get_smax()const{return *smax;}

        dTensorBC5& fetch_L() {return *L;}
        const dTensorBC5& get_L()const{return *L;}

    public:
        // additional functionality for time stepping methods
        RKinfo rk;        


    private: // data
        // data that changes but is not necessary to define the state
        // (allocated here for convenience so DogSolve* does not have
        // to reallocate it with each call)
        //
        dTensorBC4* smax;
        dTensorBC5* L;

    protected:
        virtual void initOutputDirectory(int idx, const char* framedir);

    private: // methods
        virtual void initOutputDirectories();
        virtual void saveState();
        virtual void revertToSavedState();
        int advanceTimeStepRK(double dt);

    private: // disabled
        DogSolverCart3& operator=(const DogSolverCart3& in);
        DogSolverCart3(const DogSolverCart3& in);

    private:
        const DogStateCart3& get_state_old()const
        { return (const DogStateCart3&)DogSolver::get_state_old();}
        DogStateCart3& fetch_state_old()
        {return (DogStateCart3&)DogSolver::fetch_state_old();}

    private:
        //virtual void ConstructL(const DogState& state);
        // protected: // uncomment this line if needed
        virtual double GetCFL(double dt) const;
};

void Output(const DogSolverCart3& solver, int n);
void Restart(DogSolverCart3& solver, int nstart);

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
void ConstructL(const dTensorBC5& aux,
		const dTensorBC5& q,
		dTensorBC5& Lstar,
		dTensorBC4& smax);

// ----------------------
// user-defined callbacks
// ----------------------
//
// limiter
//
void ApplyLimiter(dTensorBC5& aux, 
		  dTensorBC5& q,
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
void AfterQinit(DogSolverCart3& solver);

///////////////////////////////////////////////////////////////////////////////
// EVENTUALLY only the DogSolverCart3 argument should remain and the rest
// should be deleted. In order to do this, however, need to make some
// changes to DogSolveRK and DogSolveSDC
///////////////////////////////////////////////////////////////////////////////
void BeforeStep  (double dt, dTensorBC5& aux, dTensorBC5& q, DogSolverCart3& solver);
void AfterStep   (double dt, dTensorBC5& aux, dTensorBC5& q, DogSolverCart3& solver);
void AfterReject (double dt, dTensorBC5& aux, dTensorBC5& q, DogSolverCart3& solver);
//////
void AfterUpdateSoln  (const dTensorBC5& aux, dTensorBC5& q, double dt, double beta);
void AfterFullTimeStep(DogSolverCart3& solver);
void ConSoln          (const dTensorBC5& aux, const dTensorBC5& q, double t);

//
// called on startup
//
void InitApp();
double Qinit_restart(int nstart, dTensorBC5& q);

#endif
