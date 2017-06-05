function [XX,YY, ppApprox, ppError, er2pp,  erinfpp, uex] = postprocess2d( ...
    u, time, mpp, xbounds, ybounds )
%POSTPROCESS.    Calculate post-processed solution, and evaluate it.
%
% This version only uses the symmetric post-processor
%
% Input:
% ------
%
%   u(1:mx_old, 1:my_old, 1:kmax, 1:meqn ) : Vector of unknowns.  
%
%               No ghost cells are added to this routine. Here, 
%               kmax = # of polynomials.  It is assumed that these have
%               already been translated into (non-orthonormal) Gauss-Legendre
%               polynomials.
%
%   time  - time for solution.  This is used to construct exact solution.
%
%   mpp   - Polynomial degree.  This routine supports mpp \in { 1, 2, ... 4 }.
%
%   xbounds(1:2) - size of computational domain.
%   ybounds(1:2) - size of computational domain.
%
% Returns:
% --------
%
%   XX,YY,      : grid points after adding in 6^2 Gauss quadrature points in
%                 each element.
%
%  ppApprox     : Post-processed approximation evaluated at each grid point.
%
%  ppError      : Difference between exact and approximate at each grid point.
%
%  er2pp        : L2-norm of the post-processed solution.
%
% erinfpp       : L-inf norm of post-processed solution.
% 
% uex           : Exact solution.  Can be plotted with surf(XX,YY,uex).
%
% See also: C1to6k, get_u2d.


NUMGAUSSPTS = 6;

% Domain and problem size
[mx_old, my_old, kmax, meqn ] = size(u);

mx = NUMGAUSSPTS*mx_old;
my = NUMGAUSSPTS*my_old;

% Set up "meshinfo".  
xleft   = xbounds(1); xright  = xbounds(2);
yleft   = ybounds(1); yright  = ybounds(2);

% dx, and dy only make sense with respect to the uniform grid
dx      = (xright-xleft)/mx_old;
dy      = (yright-yleft)/my_old;

% Vector of cell centers for each element.  
x = linspace( xleft+0.5*dx, xright-0.5*dx, mx_old )';
y = linspace( yleft+0.5*dy, yright-0.5*dy, my_old )';

% Grid points (same format as meshgrid)
XX     = zeros(mx, my );
YY     = zeros(mx, my );

% Gauss-legendre quadrature points and weights
[wg, xg] = gleg_points();

% Approximation at each quadrature point and exact solution
ppApprox = zeros(mx, my, meqn);
ppError  = zeros(mx, my, meqn);
uex      = zeros(mx, my, meqn);

% Post-processed Errors 
er2pp   = zeros( 1, meqn )';
erinfpp = zeros( 1, meqn )';

% Stencil for computing errors
% shift function (Stencil for postprocesser is 
% ( x_{i\pm kwide}, y_{i\pm kwide} )
kwide = int16((3*mpp+1)/2);

% ---------------------------------------------------------------------------- %
% DEFINE THE COEFFEICIENTS used for constructing post-processed solution.
% ---------------------------------------------------------------------------- %
ccmod = C1to6k( mpp );

% ---------------------------------------------------------------------------- %
% Main loop over entire domain
% ---------------------------------------------------------------------------- %
for j=1:my_old
for i=1:mx_old

    % Loop over each Gaussian quadrature point
    for m1=1:NUMGAUSSPTS
    for m2=1:NUMGAUSSPTS

        % Index into the large grid
        ii = m1 + NUMGAUSSPTS*(i-1);
        jj = m2 + NUMGAUSSPTS*(j-1);

        % Grid location for this (fine) point
        xr = x(i)+xg(m1)*dx*0.5; 
        yr = y(j)+xg(m2)*dy*0.5;
        XX(ii,jj) = xr;
        YY(ii,jj) = yr;

        % Exact solution (evaluate the [periodic] exact solution)
        uex(ii,jj, : ) = q_init2d( xr-time, yr-time );

        % Symmetric Post-processed Approximation
        for me=1:meqn

            ppuapprox = 0.;

            k = 1;
            for degp1 = 1:(mpp+1)
            for k1  = (degp1):-1:1
                k2  = (degp1)-k1+1;

                for kgamx = -kwide:kwide
                for kgamy = -kwide:kwide
                    
                    % Hard coded boundary conditions enforced here
                    imod = 1 + mod( i+kgamx-1, mx_old );
                    jmod = 1 + mod( j+kgamy-1, my_old );
                    ppuapprox = ppuapprox +  ...
                        ccmod(kgamx+kwide+1,k1,m1)*ccmod(kgamy+kwide+1,k2,m2)*u(imod, jmod, k, me);
                end
                end
                k = k+1;
            end
            end
            ppApprox(ii,jj,me) = ppuapprox;

            % Errors (L2-error is an integral over entire domain)
            er2pp(me) = er2pp(me) + (uex(ii,jj,me) - ppuapprox)^2*(dx*dy*wg(m1)*wg(m2));
        end

    end
    end

end
end

% Pointwise (L-inf) error
ppError = uex - ppApprox;
erinfpp = max( max( abs(ppError) ) );

% L2 error
er2pp   = sqrt(er2pp*0.25);

end
