#ifndef DogStateCart4_h
#define DogStateCart4_h
#include "DogState.h"

class dTensorBC4;
class dTensorBC5;
class dTensorBC6;
class DogSolverCart4;

struct DogStateCart4: public DogStateTB
{

    private:
        dTensorBC6* q;
        dTensorBC6* aux;

    private:
        void init(int mx, int my, int mz, int mw);
        DogStateCart4& operator=(const DogStateCart4& in); // disabled
        DogSolverCart4& fetch_solver(){return (DogSolverCart4&) DogState::get_solver();}

    protected:
        virtual const char* get_classname(){return (const char*) "DogStateCart4";}
        DogStateCart4(const DogStateCart4&in, CopyMode::Enum copyMode=CopyMode::DEEP);

    public:

        // Default constructor:
        DogStateCart4(): q(NULL), aux(NULL) {};

        // Destructor:
        ~DogStateCart4();

        virtual DogState* clone(CopyMode::Enum copyMode)const
        { return new DogStateCart4(*this, copyMode); }
        virtual void init();
        virtual void copyfrom(const DogStateCart4& in);
        virtual bool check_valid_state() const;

// Reverting to top level debug_check_condition().
//     It looks like this was used once to check the posititivy of the
//     function, which I beleive is now being checked elsewhere.  
//     (-DS, 3/27/2013).
//      virtual bool debug_check_condition() const;

        virtual void init_coarser_state(const DogStateCart4& in, int mx, int my, int mz, int mw);

        // write restart frame and movie frames
        virtual void write_output(int noutput) const;
        virtual void write_frame(int nframe, const char* outputdir) const;
        // read (restart) frame
        virtual void read_frame(int nframe);

        // accessors (fetch means write-access, get means read-access)
        //
        const dTensorBC6& get_q()   const{return *q;  }
        const dTensorBC6& get_aux() const{return *aux;}
        dTensorBC6& fetch_q()            {return *q;  }
        dTensorBC6& fetch_aux()          {return *aux;}
        const DogStateCart4& get_solver()const
        {return (DogStateCart4&) DogState::get_solver();}

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
        const dTensorBC6& q, double t, int nframe);
double ReadStateArray(const char* dir, const char* varname,
        dTensorBC6& q, int nstart);

#endif
