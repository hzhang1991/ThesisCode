// Coefficients used for solving Advection equation using Lax-Wendroff. time
// stepping.
//
// When written in the form 
//
//      q_t = -u(t) q_x
//
// we have ( q(t+dt) - q(t) ) / dx = - F_x
//
// with F := F1 q + F2 q_x + F3 q_{xx} + \cdots
//
// and the coefficients are
//
// F1 = u + dt / 2 * u_t + 1/6 dt^2 u_tt
// F2 = -dt/2 * ( u^2 + dt * u * u_t )
// F3 = 1/6 * dt^2 * u^3
//
// for the third order case
struct LxW_Coeffs{

//	double dx;
	double F1;
	double F2;
	double F3;
	double F4;

	// advection speeds and the derivatives of the advection speeds
	double u;
	double ut;
	double utt;
	double uttt;
	double utttt;

};

extern LxW_Coeffs LxWcoeffs;  
