function [XX,YY, Approx, Error, er2, erinf] = ...
    construct_pre_post_process_errors2d(u, time, mpp, xbounds, ybounds )
%CONSTRUCT_PRE_POST_PRCOESS_ERRORS.  Evaluate the solution and compute errors.

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

% Gauss-legendre polynomial at each quadrature point
[wg, xg] = gleg_points();
pl       = evaluate_phi2d(mpp, xg);

% Approximation at each quadrature point and exact solution
Approx = zeros(mx, my, meqn);
uex    = zeros(mx, my, meqn);

% Errors 
er2   = zeros( 1, meqn )';
erinf = zeros( 1, meqn )';

% ---------------------------------------------------------------------------- %
% Main loop over entire domain
% ---------------------------------------------------------------------------- %
for j=1:my_old
for i=1:mx_old

    % Look at pointwise errors over grid cell (x_i, y_j)
    %
    % Note that we save a list of points, using the index m1,m2 for the
    % current cell.
    mp   = 1;
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

        % Approximation (evaluate the basis functions)
        uval = zeros(meqn);
        for me = 1:meqn
        for k  = 1:kmax
            uval(me) = uval(me) + u(i,j,k,me)*pl(k,mp);
        end
        end
        Approx(ii,jj,me) = uval(me);

        % Errors (L2-error is an integral over entire domain)
        for me =1:meqn
            er2(me) = er2(me) + 0.25*(uex(ii,jj,me) - uval(me))^2*dx*dy*wg(m1)*wg(m2);
        end
        mp = mp+1;
    end
    end

end
end

% Pointwise error for plotting.
Error = uex-Approx;

% L-inf error and L2-error
er2   = sqrt(er2);
erinf = max( max( abs(Error) ) );

end
