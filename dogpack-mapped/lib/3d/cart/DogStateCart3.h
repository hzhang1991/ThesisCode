#ifndef DogStateCart3_h
#define DogStateCart3_h
#include "DogState.h"

class dTensorBC4;
class dTensorBC5;
class DogSolverCart3;

struct DogStateCart3: public DogStateTB
{

    private:
        dTensorBC5* q;
        dTensorBC5* aux;

    private:
        void init(int mx, int my, int mz);
        DogStateCart3& operator=(const DogStateCart3& in); // disabled
        DogSolverCart3& fetch_solver(){return (DogSolverCart3&) DogState::get_solver();}

    protected:
        virtual const char* get_classname(){return (const char*) "DogStateCart3";}
        DogStateCart3(const DogStateCart3&in, CopyMode::Enum copyMode=CopyMode::DEEP);

    public:

        // Default constructor:
        DogStateCart3(): q(NULL), aux(NULL) {};

        // Destructor:
        ~DogStateCart3();

        virtual DogState* clone(CopyMode::Enum copyMode)const
        { return new DogStateCart3(*this, copyMode); }
        virtual void init();
        virtual void copyfrom(const DogStateCart3& in);
        virtual bool check_valid_state() const;

// Reverting to top level debug_check_condition().
//     It looks like this was used once to check the posititivy of the
//     function, which I beleive is now being checked elsewhere.  
//     (-DS, 3/27/2013).
//      virtual bool debug_check_condition() const;

        virtual void init_coarser_state(const DogStateCart3& in, int mx, int my, int mz);

        // write restart frame and movie frames
        virtual void write_output(int noutput) const;
        virtual void write_frame(int nframe, const char* outputdir) const;
        // read (restart) frame
        virtual void read_frame(int nframe);

        // accessors (fetch means write-access, get means read-access)
        //
        const dTensorBC5& get_q()   const{return *q;}
        const dTensorBC5& get_aux() const{return *aux;}
        dTensorBC5& fetch_q() {return *q;}
        dTensorBC5& fetch_aux() {return *aux;}
        const DogStateCart3& get_solver()const
        {return (DogStateCart3&) DogState::get_solver();}

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
