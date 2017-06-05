function LegVals = GetHybr4Legendre(kmax4d, z4d)
%GETHYBR4LEGENDRE   Sample the Hybrid 4D Legendre polynomials
%
% Usage: LegVals = GETHYBR4LEGENDRE(kmax, z4d );
%
% Sample 4D Hybrid Legendre polynomials that are on a cross-product of a
% structured and an unstructured grid.
%
% 4D-coordinates are given by a tensor product of unstructured and structured
% basis functions.
%
% The unstructured grid is on (x,y).  Canonical variables are mapped
% from:
%           (x,y) <-> (xi,eta);  
%
% The structured grid is on (vx,vy). Canonical variables are mapped
% from:
%
%         (vx,vy) <-> (mu,tau);
%
% Parameters:
%
%    kmax              - number of polynomials being sampled
%    z4d(1:mpts, 4)    - list of points that are sampled, described by:
%        z(:,1) = xi; z(:,2) = eta; z(:,3) = mu; z(:,4) = tau.
%
% Returns:
%
%    LegVals(1:kmax, 1:numpts) - The legendre polynomials evaluated 
%                                at every point indicated in the set of 
%                                points listed by z.
%
% See also: GetCart2Legendre, GetUnst2Legendre.

  switch( kmax4d )
    case 15
      spatial_order = 3;
    case 5
      spatial_order = 2;
    case 1
      spatial_order = 1;
    otherwise
      disp('Unsupported spatial order provided');
      exit(1);
  end

  % Find out how many points we want to sample:
  [mpts, ndim] = size( z4d );
  assert( ndim == 4 );

  kmax2d  = get_kmax( spatial_order, 2 );

  spts_unst    = [z4d(:,1) z4d(:,2)];
  spts_cart    = [z4d(:,3) z4d(:,4)];

  if( spatial_order > 2 )
      LegVals_Unst = GetUnst2Legendre( kmax2d, spts_unst );
      LegVals_Cart = GetCart2Legendre( kmax2d, spts_cart );
  end

  % Canonical variables
  xi  = z4d( :, 1 );
  eta = z4d( :, 2 );
  mu  = z4d( :, 3 );
  tau = z4d( :, 4 );

  % Placeholder for Legendre values
  LegVals = ones( kmax4d, mpts );

  sq2 = sqrt(2.0);  sq3 = sqrt(3.0);
  if( spatial_order > 1 )

    LegVals( 5, : ) = sq3*tau(:);
    LegVals( 4, : ) = sq3* mu(:);
    LegVals( 3, : ) = sq2*sq3*( xi(:) + 2*eta(:) );
    LegVals( 2, : ) = 3.0*sq2*xi(:);

  end

  if( spatial_order > 2 )

    % 'Cartesian' terms:
    LegVals( 15, :) = LegVals_Cart(6, :);
    LegVals( 14, :) = LegVals_Cart(5, :);
    LegVals( 13, :) = LegVals_Cart(4, :);

    % cross (cartesian + unstructured) terms:
    LegVals( 12 , :) = LegVals_Unst(3, :).*LegVals_Cart(3, :);
    LegVals( 11 , :) = LegVals_Unst(3, :).*LegVals_Cart(2, :);
    LegVals( 10 , :) = LegVals_Unst(2, :).*LegVals_Cart(3, :);
    LegVals(  9 , :) = LegVals_Unst(2, :).*LegVals_Cart(2, :);

    % 'unstructured' terms:
    LegVals( 8, :) = LegVals_Unst(6, :);
    LegVals( 7, :) = LegVals_Unst(5, :);
    LegVals( 6, :) = LegVals_Unst(4, :);

  end

end
