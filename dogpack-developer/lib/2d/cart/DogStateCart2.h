#ifndef DogStateCart2_h
#define DogStateCart2_h
#include "DogState.h"

class dTensorBC3;
class dTensorBC4;
class DogSolverCart2;

class DogStateCart2: public DogStateTB
{

    private:
        dTensorBC4* q;
        dTensorBC4* aux;

    private:
        void init(int mx, int my);
        DogStateCart2& operator=(const DogStateCart2& in); // disabled
        DogSolverCart2& fetch_solver(){return (DogSolverCart2&) DogState::get_solver();}

    protected:
        virtual const char* get_classname(){return (const char*) "DogStateCart2";}
        DogStateCart2(const DogStateCart2&in, CopyMode::Enum copyMode=CopyMode::DEEP);

    public:

        // Default constructor:
        DogStateCart2(): q(NULL), aux(NULL) {};

        // Destructor:
        ~DogStateCart2();

        virtual DogState* clone(CopyMode::Enum copyMode)const
        { return new DogStateCart2(*this, copyMode); }
        virtual void init();
        virtual void copyfrom(const DogStateCart2& in);
        virtual bool check_valid_state() const;

// Reverting to top level debug_check_condition().
//     It looks like this was used once to check the posititivy of the
//     function, which I beleive is now being checked elsewhere.  
//     (-DS, 3/27/2013).
//      virtual bool debug_check_condition() const;

        virtual void init_coarser_state(const DogStateCart2& in, int mx, int my);

        // write restart frame and movie frames
        virtual void write_output(int noutput) const;
        virtual void write_frame(int nframe, const char* outputdir) const;
        // read (restart) frame
        virtual void read_frame(int nframe);

        // accessors (fetch means write-access, get means read-access)
        //
        const dTensorBC4& get_q()   const{return *q;}
        const dTensorBC4& get_aux() const{return *aux;}
        dTensorBC4& fetch_q() {return *q;}
        dTensorBC4& fetch_aux() {return *aux;}
        const DogStateCart2& get_solver()const
        {return (DogStateCart2&) DogState::get_solver();}

        // -- user "callbacks" -- //
        //
        // class methods with the virtual keyword have subclasses that can override
        // these methods.  In other words, EVERY subclass of this class will have
        // these functions, and can override them.

        // This function projects initial conditions onto aux and q.
        virtual void InitState();

        virtual void AfterInitState();
        virtual void BeforeStage(double dt);
        virtual void SetBndValues();
        virtual void ConstructL();
        virtual void ApplyLimiter();
        virtual void AfterStage(double dt);
        virtual void AfterReject(double dt);
        virtual void AfterFullTimeStep(double dt);
        virtual void ReportAfterStep() const;
};

// low-level methods to write and read state arrays
void WriteStateArray(const char* framedir, const char* varname,
        const dTensorBC4& q, double t, int nframe);
double ReadStateArray(const char* dir, const char* varname,
        dTensorBC4& q, int nstart);

#endif
