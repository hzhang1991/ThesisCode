function [u, time] = get_u2d( outputdir, frameNum, Elmts, meqn, morder )
%GET_U2D.  Function to read u from a file.
%
% This function was written to be able to interface with DoGPack, where the
% polynomials are rescaled to be orthonormal.  Additionally, DoGPack saves the
% final time of the solution.
%
% Parameters
% ----------
%
%           fname - name of file to extract data from
%
%           Elmts - number of grid elements for that solution.  (This is used
%           to reshape the vector of unkowns.)
%               mx = Elmts(1), my = Elmts(2).
%
%           morder - order of accuracy.  (morder == mpp+1)
%
% Returns
% -------
%
%           u - vector of unkowns. u = u( 1:mx, 1:my, 1:kmaxbig, 1:meqn )
%               NOTE: this is a vector of (nonorthonormal) Legendre
%               coefficients!  They have been converted from what DoGPack
%               first constructs.
%
%               kmaxbig = (mpp+1)^2.  DoGPack does not save extra coefficients,
%               here we compute the full tensor product because it is easier
%               to code!  Extra coefficients are set to zero.
%
%       time  - time of solution.  DoGPack saves the time for each output, so
%               we may as well return that as well.
%
% See also: main_driver, and 
%   $DOGPACK/lib/2d/cart/Legendre2d.cpp:SetLegendreAtPoints
% for how the coefficients are ordered.

    frameNum = int32( frameNum );
    mx = Elmts(1); my = Elmts(2);

    % Read in the DoGPack coefficients.
    % NOTE: size(q_data) = [mx,my,kmax, 1:meqn]
    kmax = get_kmax(morder, 2);
    [q_data,time] = read_state2_cart( 1, outputdir, frameNum, 'q', mx, my, meqn, kmax, 1:meqn);
    u = q_data;

    % Reorder the third-order terms, because these were not printed in a
    % 'tensor product' order.
    % 
    % See: SetLegendrePolys in $DOGPACK/lib/2d/cart/L2Project.cpp for a list
    % of the basis functions we use.
    %
    % The new polynomials will be order in the following manner, before
    % orthogonalization:
    %
    %    { 1, xi, eta, xi^2, xi*eta, eta^2, xi^3, xi^2*eta, xi*eta^2, eta^3, ... }
    %
    switch( morder )

        case 1
        P = [1];

        case 2
        P = [1,2,3];

        case 3
        P = [1,2,3,5,4,6];

        case 4
        P = [1,2,3,5,4,6,9,7,8,10];

        case 5
        P = [1,2,3,5,4,6,9,7,8,10,14,11,13,12,15];

    end
    q_ordered = q_data(:,:,P,:);

    % We will zero out every higher-order term in the Galerkin expansion.
    u = zeros( mx, my, morder^2, meqn );
    u(:,:,1:kmax,1:meqn) = q_ordered(:,:,:,:);

    % Rescale each polynomial in u
    for me=1:meqn
    for j=1:my
    for i=1:mx
        k = 1;
        for degp1 = 1:morder
        for k1  = (degp1):-1:1
            k2  = (degp1)-k1+1;
            u(i,j,k,me) = u(i,j,k,me) * sqrt( 2.0*k1-1.0 ) * sqrt(2.0*k2-1.0 );
            k = k+1;
        end
        end
    end
    end
    end

end
