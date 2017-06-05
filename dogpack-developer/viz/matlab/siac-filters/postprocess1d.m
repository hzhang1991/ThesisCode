% ---------------------------------------------------------------------------- %
function [Ypp, ppApprox, ppError, er2pp,  erinfpp] = postprocess1d( u, time, mpp, xbounds, Elmts )
%POSTPROCESS.    Calculate post-processed solution, and evaluate it.
%
% This version only uses the symmetric post-processor
%
% Input:
% ------
%
%   u(mpp+1,  Elmts+2 ) - vector of unknowns.  u(:,1) and u(:,Elmts+2) are ghost
%                         cells
%
%   time  - time for solution.  This is used to construct exact solution.
%
%   mpp   - Polynomial degree.  This routine supports mpp \in { 1, 2, ... 4 }.
%
%   xbounds(1:2) - size of computational domain.
%
%   Elmts - Number of grid elements.  The interior of the domain is i=2:Elmts+1.
%                                     i=1, and i=Elmts+2 are ghost elements.
%
% Returns:
% --------
%
%   er2pp           - l2 error
%   erinfpp         - l-inf error
%
% See also: lineartransport or lineartransport_nopp for the main drivers. and 
%           C1to6k[1-4] for defining coefficients (used for B-splines?).
% ---------------------------------------------------------------------------- %

% Set up the meshinfo.
xleft = xbounds(1); xright = xbounds(2);
dx = (xright-xleft)/Elmts;

% Vector of cell centers for each element.  This adds a single ghost cell to
% the left and right of the domain.
x = linspace( xleft-0.5*dx, xright+0.5*dx, Elmts+2 )';  % Same as above for loop

morder  = mpp+1;

% shift function
kwide = int16((3*mpp+1)/2);

% ---------------------------------------------------------------------------- %
%DEFINE THE COEFFEICIENTS used for constructing post-processed solution.
% ---------------------------------------------------------------------------- %
ccmod = C1to6k( mpp );

% ---------------------------------------------------------------------------- %
% Compute errors when compared with exact solution (sampled at GL-points)
% ---------------------------------------------------------------------------- %
er2pp    = 0; erinfpp  = 0;

% Each of these will eventually have size ( (Elmts+1)*6, 1 )
Ypp      = [];      % Mesh - sampled at each quadrature point
Exactpp  = [];      % Exact solution - sampled at each quadrature point
ppApprox = [];      % Symmetric, post-processed solution sampled at each quadrature point

% Gaussian quadrature weights and points, (using exactly 6 points)
[wg, xg] = gleg_points();

for i=2:Elmts+1  % modified to loop over interior points only (-DS)

    % Loop over each Gaussian quadrature point
    for j=1:6
        xr=x(i)+xg(j)*dx*0.5;
        xrt=xr-time;
        
        % Symmetric Post-processed Approximation
        ppuapprox=0;
        for k=1:morder
            for m=-kwide:kwide

                % Hard-coded periodic boundary conditions
                if ((i+m)<=0)
                    ppuapprox = ppuapprox + ccmod(m+kwide+1,k,j)* ...
                            u(k,Elmts+(i+m));
                elseif ((i+m)>=(Elmts+2))
                    ppuapprox = ppuapprox + ccmod(m+kwide+1,k,j)* ...
                            u(k,2+i+m-(Elmts+2));
                else
                % "Interior" elements
                    ppuapprox = ppuapprox + ccmod(m+kwide+1,k,j)*u(k,i+m);
                end    
            end
        end
        ppApprox=[ppApprox;ppuapprox];

        % Exact Solution
        exact   = q_init1d( xrt );
        Exactpp = [Exactpp;exact];
        Ypp     = [Ypp;xr];

        % Running total for l2-error
        er2pp=er2pp+(exact-ppuapprox)^2*dx*wg(j);

    end
end

% Pointwise (L-inf) error
ppError = abs(Exactpp-ppApprox);
erinfpp = max(ppError);

% L2 error
er2pp   = sqrt(er2pp/2);

end
