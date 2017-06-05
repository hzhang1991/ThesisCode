#include "DogState.h"

// -------------------------------------------------------------------------- //
// Currently, the only State class that overrides this is DogStateCart2.
// There, it was used to see if a positivity violation occurred.
//
// Returns:
//
//      True  : if in the all clear to keep moving forward.
//
//      False : if the solution is corrupt (e.g. violated positivity).
// -------------------------------------------------------------------------- //
bool DogState::debug_check_condition() const
{ return true; }

void DogState::InitState(){time=0;}

// default implementations of virtual methods
int DogState::advanceSplitTimeStep(double dt)
{
    eprintf("splitting is turned on but no split operator is defined.\n"
            "\tApplications that use splitting must override this function\n"
            "\tor the functions that by default call it.");
    return 0;
}
