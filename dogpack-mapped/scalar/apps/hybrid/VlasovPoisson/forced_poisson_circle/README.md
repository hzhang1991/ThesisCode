Forced Problem
==============

This problem has not yet be written.  The idea is to test the Poisson solver,
one dimension at a time.

Poisson's Equation
------------------

Take the electron density to be:

    rho = 1-x^2-y^2

and then what follows is the solution to Poisosn's equation is given by:

    phi = x^4/12 + y^4/12 -x^2/4 - y^2/4

Both components of the electric field are then:

    E1 = -phi_x = x/2 - x^3/3,
    E2 = -phi_y = y/2 - y^3/3.

Source Term
-----------

We use the method of manufactured solutions:

    f_t + vx f_x + vy f_y - E1 f_vx - E2 f_vy = psi,

and choose an exact solution f which provides these values.

Example 1
---------

Select the exact solution to be:

    f = ( 3/4 + (1/4)*cos(2*Pi*t) ) * (1-x^2-y^2) * exp( -vx^2 ) / sqrt(Pi);

This doesn't depend on vy, so E2 doesn't show up in the source term.  The
exact density is given by,

    rho = ( 3/4 + (1/4)*cos(2*Pi*t) ) * (1-x^2-y^2),

and the exact electric field is given by:

    E1 = -phi_x = (3/4 + 1/4 cos(2*pi*t) ) * (-x/2 + x^3 / 3 ),
    E2 = -phi_y = (3/4 + 1/4 cos(2*pi*t) ) * (-y/2 + y^3 / 3 ).

