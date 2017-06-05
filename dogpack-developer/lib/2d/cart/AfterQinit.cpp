#include<iostream>
using namespace std;

class DogSolverCart2;

// Routine that is called *after* the initial conditions are projected onto the
// basis functions.  This routine is useful, for example, for pre-processing of
// initial data, saving information about the initial conditions (e.g. initial
// density or momentum) or a variety of other options.
//
// This default routine does nothing.  Modify your personal Makefile if you wish
// to use this function call.
//
// Input:
// ------
//
//      DogSolverCart2 solver      - Large wrapper for 2D Cartesian solver.
//                                   This gives access to current state, q and
//                                   other items.
//
// Returns:
// --------
//
//      nothing - unless you write this routine!
//
// See also: ...
void AfterQinit(DogSolverCart2& solver)
{

/*
 * Here is an example of accessing the usual stuff.  Some extra items need to be
 * added to the above to uncomment this and make this compile.
 *

 *  DogStateCart2* dogStateCart2 = &solver.fetch_state();
 *  dTensorBC4& qnew = solver.fetch_state().fetch_q();
 *  dTensorBC4& aux  = solver.fetch_state().fetch_aux();

 *  const int mx     = qnew.getsize(1);
 *  const int my     = qnew.getsize(2);
 *  const int maux   = aux.getsize(3);
 *  const int meqn   = qnew.getsize(3);
 *  const int kmax   = qnew.getsize(4);
 *  const int mbc    = qnew.getmbc();

 */

}
