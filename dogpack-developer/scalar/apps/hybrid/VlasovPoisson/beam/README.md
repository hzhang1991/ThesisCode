Beam problem from:

Besse, Segre and Sonnendrucker, "Semi‐Lagrangian Schemes for the Two‐Dimensional Vlasov‐Poisson System on
Unstructured Meshes", Transport Theory and Statistical Physics (2005).

See:

    http://dx.doi.org/10.1080/00411450500274592

The initial conditions are of the form:

    f0 = n0 / (2 pi vth^2 pi a^2 ) exp( -(vx^2+vy^2)/( 2 vth^2 ) ) exp ( -(x^2+y^2)/( 2 R^2 ) )

The parameters we use are:

    a   = 1
    vth = 1
    R^2 = 1/4

    eta    = 1/4
    omega0 = omega / eta
    n0     = 2*(1-eta^2) * omega0^2

TODO - the paper doesn't say what to set omega equal to.  I will set omega0 = 1,
and that implies omega = 1/4.
