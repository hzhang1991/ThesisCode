#include "DogStateUnst2.h"

DogStateUnst2 dogStateUnst2;

DogStateUnst2::DogStateUnst2()
{
    time=0.0;
    dt=0.0;
    initial_dt=0.0;
}

DogStateUnst2::DogStateUnst2(DogStateUnst2& old_dogStateUnst2)
{
    time       = old_dogStateUnst2.time;
    dt         = old_dogStateUnst2.dt;
    initial_dt = old_dogStateUnst2.initial_dt;
}

// the 'do nothing' destructor
DogStateUnst2::~DogStateUnst2()
{
}

void DogStateUnst2::init()
{
    time       = 0.0;
    dt         = 0.0;
    initial_dt = 0.0;
}

