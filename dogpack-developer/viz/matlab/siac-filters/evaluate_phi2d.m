function pl = evaluate_phi2d(mpp, xg)
%EVALUATE_PHI2d.  Returns the basis functions evaluated at each quadrature point.
%
% Input:
% ------
%
%     mpp : degree of polynomial.  mpp \in {0, 1, ..., whatever }.
%     xg  : Quadrature points to evaluate the polynomial
%
% Returns:
% --------
%
%     pl  : Tensor product of basis functions evaluated at each point in xg.
%     size(pl) = [ (mpp+1)^2, 6^2 ]
%
% See also: gleg_points.m

    NUMGAUSSPTS = 6;

    pl   = zeros( (mpp+1)^2, NUMGAUSSPTS*NUMGAUSSPTS );
    pl1d = evaluate_phi1d( mpp, xg );

    % NOTE: order in how polynomials are evaluated!!!
    % 
    % In DoGPack, we group them by degrees.  That is,
    %
    % deg0 = {1}
    % deg1 = {xi, eta}
    % deg2 = {xi^2, xi*eta, eta^2 }  (After swapping two values ... )
    % deg = {xi^3, xi^2*eta, xi*eta^2, eta^3 }
    %
    % and so on.
    k = 1;
    for degp1 = 1:(mpp+1)
    for k1  = (degp1):-1:1
        k2  = (degp1)-k1+1;
        mp = 1;
        for m1=1:6
        for m2=1:6
            pl(k,mp) = pl1d(k1,m1)*pl1d(k2,m2);
            mp = mp+1;
        end
        end
        k = k+1;
    end
    end

end
