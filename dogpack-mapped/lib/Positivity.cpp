#include <float.h>          // for DBL_MAX
#include "debug.h"
#include "assert.h"
#include "DogParams.h"
#include "Positivity.h"

// Notes on positivity
//
// Idea: let the code run positivity is violated
//
// Positivity modes:
// + off
//   + print error and quit if cell average or riemann
//     points violate positivity
// + unconstrained
//   + fall back on stability CFL condition
//   + throw exception if cell average or riemann
//     points violate positivity
// + corrected:
//   + enforce machine positivity within cell
// + guaranteed:
//   + enforce machine-positivity within cell and
//   + use positivity-preserving time step
//
// The exception handler should:
// + increase the positivity stringency
//   + increase mode stringency if possible
//   + else increase positivity CFL stringency
// + revert the dogState
// + compute a positivity-preserving time-step
//
// In case the exception handler is not called,
// the positivity checker should
// + diminish the positivity stringency if appropriate
//   + if positivity limiters are consistently
//     not being triggered then try turning them off
//   + if cell average is never being reset then
//     relax toward less stringent CFL.
//
// I see two approaches to dealing with positivity violations:
// (1) bail out (throw an exception) when a problem arises and
// (2) systematically check positivity.
//
// The first case requires thread-safe exception handling. You
// have to catch the exception within the thread, then tell
// the other threads to stop and bail out, then probably throw
// another exception, and then catch it at the high level and
// deal with it. If you are compiling with exceptions then
// once this framework is in place the user merely must throw
// exceptions when appropriate and provide a way to handle the
// exception (e.g. by increasing positivity stringency). If you
// are not compiling with exceptions then exception-handling must
// be done "C-style" by returning error codes, which can require
// systematic code revision.
//
// The second case avoids bailing out in the middle of a
// parallelized loop through the mesh. This allows for greater
// ease of porting to parallel systems (e.g. GPU) where bailing
// out of a thread can be problematic, but portends programming
// and/or computational expense. Specifically, it requires
// (a) implementing a separate positivity-checking routine or
// (b) implementing a mechanism to register that a positivity
//     violation has occurred.
// Case (a) incurs additional computational expense.
// Case (b) can be difficult to program because of the
//   need to communicate positivity violation information
//   back up the stack.

PositivityState::PositivityState():
    cflFactor(1.),
    stringencyLevel(POSITIVITY_POINTS),
    minval(DBL_MAX),
    stepsWithoutViolation(0),
    repeatStageFlag(false),
    repeatStepFlag(false)
{}

void PositivityState::write_state(void* in)const
{
    FILE* file = (FILE*) in;
    fprintf(file,"cflFactor = %24.16e\n", cflFactor);
}

void PositivityState::read_state(void* in)
{
    FILE* file = (FILE*) in;
    fscanf(file,"cflFactor = %lf\n", &cflFactor);
}

//bool PositivityState::using_split_source() const
//{
//  if(plasmaParams.get_source_type()==SourceType::SPLIT)
//    return true;
//  if(stringencyLevel==IMPLICIT_SOURCE)
//    return true;
//  return false;
//}

// This way of handling positivity violations is
// not continuous in time.  One would have to
// transition gradually from one level to another.
// stringencyLevel should be changed to double
// (integer values could serve as points where
// complete transition is made).
//
void PositivityState::increaseStringency()
{
    if(stringencyLevel==AVERAGE)
    {
        dprintf("increasing stringency level to POSITIVITY_POINTS");
        stringencyLevel = POSITIVITY_POINTS;
        repeatStageFlag = true;
    }
    //else if(stringencyLevel==POSITIVITY_POINTS)
    //{
    //  dprintf("increasing stringency level to IMPLICIT_SOURCE");
    //  stringencyLevel = IMPLICIT_SOURCE;
    //  repeatStepFlag = true;
    //}
    //else if(stringencyLevel==IMPLICIT_SOURCE)
    else if(stringencyLevel==POSITIVITY_POINTS)
    {
        // make this number configurable?
        // could rescale based on comparison
        // of minval with minval computed with
        // larger time step.
        //
        // want to modify dt by a factor that will
        // cause minval to be zero.
        //
        // time step must be repeated if step length is changed.
        repeatStepFlag = true;
#if 0
        if(suggested_dt_changeFactor < .8)
            dt_changeFactor = .8;
        else if(suggested_dt_changeFactor > .95)
            dt_changeFactor = .95;
        else
            dt_changeFactor = suggested_dt_changeFactor;
#endif
        //double dt_changeFactor = .95;
        //cflFactor *= dt_changeFactor;
        cflFactor -= .0625;
        //dt *= dt_changeFactor;
        //dprint(dt);
        dprintf1("decreased cflFactor to %f",cflFactor);
        // if cflFactor gets too small something must have gone wrong.
        // There should exist a positivity-guaranteeing CFL number
        // that is independent of the solution (assuming that there
        // is a positivity-guaranteeing CFL number that is independent
        // of the solution for the HLLE method, which maybe isn't quite
        // true because we don't have a perfect way to obtain an upper
        // bound on physical wave speeds for the Riemann problem between
        // two cell states (Einfeldt's prescription is a not a guarantee
        // for all strictly hyperbolic systems).
        //assert_gt(cflFactor, .08);
        assert_gt(cflFactor, 0.);
    }
    else
    {
        invalid_value_error(stringencyLevel);
    }
}

void PositivityState::update(double minval_in, double dt_in)
{
    double d_minval = minval_in - minval;
    double d_dt = dt_in - dt;

    if(d_minval!=0)
    {
        double ratio = d_dt/d_minval;
        // assumes an approximately linear relationship between dt and minval
        double suggested_dt = dt-d_dt/d_minval*minval;
        double suggested_cflChangeFactor = (suggested_dt/dt)*.99;
        if(minval_in < 0.)
            dprint1(suggested_cflChangeFactor);
    }

    set_minval(minval_in);
    set_dt(dt_in);
    if(minval_in<0.)
    {
        increaseStringency();
        stepsWithoutViolation=0;
    }
    else
    {
        stepsWithoutViolation++;
    }
    //if(stepsWithoutViolation>=1024)
    //{
    //  dprint1("resetting stringency");
    //  resetStringency();
    //}
}
