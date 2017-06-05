#include "DogStateHybrid.h"

DogStateHybrid dogStateHybrid;

DogStateHybrid::DogStateHybrid()
{
    time=0.0;
    dt=0.0;
    initial_dt=0.0;
}

DogStateHybrid::DogStateHybrid(DogStateHybrid& old_dogStateHybrid)
{
    time       = old_dogStateHybrid.time;
    dt         = old_dogStateHybrid.dt;
    initial_dt = old_dogStateHybrid.initial_dt;
}

// the 'do nothing' destructor
DogStateHybrid::~DogStateHybrid()
{
}

void DogStateHybrid::init()
{
    time       = 0.0;
    dt         = 0.0;
    initial_dt = 0.0;
}

