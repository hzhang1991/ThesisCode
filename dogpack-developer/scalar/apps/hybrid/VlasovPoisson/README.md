This folder contains all the Vlasov examples.

Boundary Conditions
===================

The library called lib/ contains the default boundary conditions.  As of now,
we currently only support the following boundary conditions for the E-field,
and density function (in phase space):

* Dirichlet on Electric Field 
* Dirichlet on density function
* periodic in vx and vy is hard coded in the semi-Lagrangian step on velocity
  space. That function is called StepAdvecCC.cpp and is located in lib/hybrid.
   

List of Examples
================

circ
----

Initial conditions: Maxwellian with uniform density.  Calls Poisson solver, and
the whole works.

test\_intitial\_conditions
------------------------

The purpose of this routine is to test the output mechanisms.
As of the writing of this Makefile, this application only prints density.
(copied from circ example).


circle\_w\_source
-----------------

Forced advection problem:

    q_t + vx q_x + vy q_y = Psi

psi is chosen so we know what the exact solution is. Grid refinement is
performed on the unstructured grid only.

force\_poisson\_circle
----------------------

Forced advection problem:

    q_t + vx q_x + vy q_y - E1 q_{vx} = Psi

In this example, E1 is computed exactly, and no Poisson solver is ever called.


