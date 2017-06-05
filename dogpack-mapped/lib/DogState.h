#ifndef _DOG_STATE_H_
#define _DOG_STATE_H_

#include "debug.h"
#include "tensors1d.h" // for CopyMode

// The descendents of the DogState class are expected
// to contain all the time-dependent state information
// of the solver, so that if a time step must be redone
// we can simply revert to the old state.
//
// To avoid casting each instance of e.g. DogStateCart2 could
// have a reference to DogSolverCart2.  Alternatively,
// we can cast under the hood.
//
class DogSolver;

struct DogState
{

    public:

        // Constructor and destructor
        DogState():time(0),solver(NULL){}
        virtual ~DogState(){}

        // user is responsible to delete output
        virtual DogState* clone(CopyMode::Enum copyMode)const=0;

        // if we made this a virtual method
        // then we would have to cast the argument
        void copyfrom(const DogState& in){time = in.time;}

        virtual bool check_valid_state() const=0;

        // a mechanism to allow the user to check an application-
        // specific condition at any point in the code when debugging
        virtual bool debug_check_condition() const;

        // accessors (fetch means write-access, get means read-access)
        //
        double    get_time() const                      { return time; }
        double& fetch_time()                            { return time; }
        void set_time(double in)                        { time = in;   }
        const DogSolver& get_solver() const             { return *solver;}

        // could make this private if DogSolver were a friend
        void set_solver(DogSolver& in){solver = &in;}

        // For debugging, to show what this class actually is.
        virtual const char* get_classname(){return (const char*) "DogState";}

    public:

        // user "callbacks"
        //
        // The dt argument in some of these callbacks is not
        // strictly necessary, since the user can access it via
        // get_solver().get_dt(), but I include in functions where
        // the user is likely to want it; I want to avoid propagating
        // solver accesses without need. -eaj
        //
        virtual void InitState();
        virtual void AfterInitState()=0;
        virtual void BeforeStage(double dt){};
        virtual void SetBndValues()=0;
        virtual void ConstructL()=0;

        // called after a step has been rejected
        virtual void AfterReject(double dt){};

        // called before limiters would be applied
        virtual bool AfterAdvanceStage(){};
        virtual void ApplyLimiter()=0;

        // called after limiters are applied
        virtual void AfterStage             (double dt){};
        virtual void AfterFullTimeStep      (double dt){};
        virtual void ReportAfterStep()const{};

        // Is this really an appropriate spot to put a time splitting
        // mechanism?  I only ask because I'm having difficulty following the
        // hierarchy of these classes.
        //
        // If I understand correctly, this class is supposed to be the 
        // "top"-level class that everyone uses to define state values.
        // Why do descendents have to inherit split time stepping
        // routines?  As far as I can tell, only the twofluid codes
        // actually use these routines, so wouldn't it be better to place it
        // elsewhere? (-DS)
        virtual int advanceSplitTimeStep    (double dt);

        int advanceFullTimeStep(double dt);

    protected:
        DogSolver& fetch_solver() { return *solver;}
        DogState(const DogState& in):time(in.time),solver(in.solver){}

    private:

        double time;
        // This is a little incestuous, but convenient to avoid
        // passing the solver around or insisting that it be a
        // singleton. Mirroring this pointer as e.g. AppStateCart2::solver
        // would additionally allow us to avoid casting in the accessors.
        DogSolver* solver;

};

// This is a simple class that stores two dTensorBase pointers: q and aux.
//
// The fecth and get routines allow one access to these values.
//
// What's the acronym TB stand for? (-DS)
class DogStateTB : public DogState
{

    public:

        // get means read access; fetch means write access
        const dTensorBase& get_q()   const{return *q;}
        const dTensorBase& get_aux() const{return *aux;}
        dTensorBase& fetch_q(){return *q;}
        dTensorBase& fetch_aux(){return *aux;}

        //bool test_aux()

        // for debugging, to show what this class actually is.
        virtual const char* get_classname(){return (const char*) "DogStateTB";}

    // "protected members are accessible within the class and its methods and
    // in its descendants"
    // http://en.cppreference.com/w/cpp/language/access
    protected: // disabled  (?? what does disabled mean ?? -DS )
        DogStateTB(const DogStateTB& in):DogState(in){}

    protected:

        // Constructor - arrays initialized to be NULL pointers.
        DogStateTB():q(0),aux(0){}

        // Methods for redirecting the pointers to q and aux to new
        // dTensorBase objects.
        void set_q(dTensorBase* in){q=in;}
        void set_aux(dTensorBase* in){aux=in;}

    private:

        // DogStateTB does not own these.  (would require making ~dTensorBase virtual.)
        //
        // [ I don't understand this comment -DS ]
        dTensorBase* q;
        dTensorBase* aux;

};

#endif
