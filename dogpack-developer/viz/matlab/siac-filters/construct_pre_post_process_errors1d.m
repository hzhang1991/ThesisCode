function [Y, Approx, Error, er2, erinf] = construct_pre_post_process_errors1d(u, time, mpp, xbounds, Elmts)
%CONSTRUCT_PRE_POST_PRCOESS_ERRORS.  Evaluate the solution and compute errors.
%
% See also: postprocess1d.m.

% Set up "meshinfo".  
xleft   = xbounds(1); xright  = xbounds(2);
dx = (xright-xleft)/Elmts;

% Vector of cell centers for each element.  
% NOTE: This adds a single ghost cell to the left and right of the domain.
x = linspace( xleft-0.5*dx, xright+0.5*dx, Elmts+2 )';

Approx = [];        % Pointwise approximation to solution
Exact  = [];        % Exact solution
Y      = [];        % Grid points
er2=0; erinf = 0;   % errors

% Gauss-legendre polynomial at each quadrature point
[wg, xg] = gleg_points();
pl = evaluate_phi1d(mpp, xg);

% ---------------------------------------------------------------------------- %
% Main loop over entire domain
% ---------------------------------------------------------------------------- %
for i=2:Elmts+1

    % Look at pointwise errors over grid cell x_i
    for j=1:6
        xr=x(i)+xg(j)*dx*0.5;
        xrt=xr-time;

        % Approximation (evaluate the basis functions)
        uapprox=0;
        for kk=1:mpp+1
            uapprox=uapprox+u(kk,i)*pl(kk,j);
        end
        Approx=[Approx;uapprox];

        % Exact solution (evaluate the [periodic] exact solution)
        exact = q_init1d( xrt );
        Exact = [Exact;exact];
        Y     = [Y;xr];
        er2   = er2+(exact-uapprox)^2*dx*wg(j);
    end
end

% Pointwise error for plotting.
Error = abs(Exact-Approx);

% L-inf error and L2-error
erinf = max(Error);
er2   = sqrt(0.5*er2);

end
