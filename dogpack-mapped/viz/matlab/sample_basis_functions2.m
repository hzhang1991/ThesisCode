function phi = sample_basis_functions2(method_order, point_type, numpts)
%SAMPLE_BASIS_FUNCTIONS2    Sample the 2D basis functions.
%
% phi = SAMPLE_BASIS_FUNCTIONS2( method_order, point_type, numpts ) samples
% the 2D basis functions at either regular intervals or Gaussian quadrature
% points.
%
% Input:
%
%   method_order = order of the method.  (See also: get_kmax)
%
%   point_type = 1:   uniform points on each element
%              = 2:   Gauss-Legendre (quadrature) points on each element
%
%   numpts - number of points to sample. ( Default: method_order )
%
% Output:
%
%   phi(1:numpts, 1:numpts, 1:kmax ) - Basis functions sampled at each point.
%
% See also: GETCART2LEGENDRE, get_kmax, sample_cart_basis_at_pts2, get_positions_1d

  if(~exist('numpts'))
    numpts=method_order;
  end

  % Find the 1D points that will be sampled:
  positions_1d = get_positions_1d(method_order, point_type, numpts);

  % Tensor product of 1D points:
  xi  = transpose(positions_1d)*ones(1,numpts);
  eta = ones(numpts,1)*positions_1d;
  phi = sample_cart_basis_at_pts2(xi, eta, method_order);

end

% point_type: 1=equispaced, 2=gaussian_quadrature
function positions_1d = get_positions_1d(method_order, point_type, numpts)
  if(point_type==1)
    % create array of
    % evaluation of Legendre basis functions at
    % method_order(times)meth1sample points in generic cell
    %
    positions_1d=-1+(2*(1:numpts)-1)/numpts;
  elseif(point_type==2)
    assert(numpts==method_order)
    if(method_order==1)
      positions_1d = 0.0;
    elseif(method_order==2)
      positions_1d = [-1.0, 1.0]/sqrt(3);
    elseif(method_order==3)
      positions_1d = [-sqrt(3), 0.0, sqrt(3)]/sqrt(5);
    elseif(method_order==4)
      p1 = sqrt(3.0+sqrt(4.8));
      p2 = sqrt(3.0-sqrt(4.8));
      positions_1d = [-p1,-p2,p2,p1]/sqrt(7);
    elseif(method_order==5)
      p1 = sqrt(5.0 + 2.0*sqrt(10)/sqrt(7));
      p2 = sqrt(5.0 - 2.0*sqrt(10)/sqrt(7));
      positions_1d = [-p1,-p2,0,p2,p1]/3.0;
    else
      error(['invalid method_order: ' num2str(method_order)]);
    end
  else
    error(['invalid point_type: ' num2str(point_type)]);
  end
end

function phi = sample_cart_basis_at_pts2(xi, eta, method_order)
%SAMPLE_CART_BASIS_AT_PTS2   Sample the 2D Cartesian basis functions.
%
% phi = SAMPLE_CART_BASIS_AT_PTS2( xi, eta, method_order ) - sample the basis
% functions at all the points (xi, eta).
%
% This function was pulled from an M-file because no-one else was calling it
% anywhere in DoGPacK.  (-DS)
%
% See also: GETCART2LEGENDRE.  These functions could be consolidated.

  % Quick error check:
  %assert( kmax==method_order*(method_order+1)/2 );
  assert(method_order<=5);

  % phi values:
  phi=zeros(size(xi,1),size(eta,2),get_kmax(method_order));

  % Quick lookup for quadrature points
  xi2=xi.*xi;
  xi3=xi.*xi2;
  xi4=xi.*xi3;
  eta2=eta.*eta;
  eta3=eta.*eta2;
  eta4=eta.*eta3;

  % Quick lookup for sqrt functions:
  sq3=sqrt(3.);
  sq5=sqrt(5.);
  sq7=sqrt(7.);

  % Each order at least retains cell averages:
  assert(method_order>=1); % 1st order
    phi(:,:,1) = 1.0;

  % Higher order basis functions:
  if(method_order>=2) % 2nd order
    phi(:,:,2) = sq3*xi;
    phi(:,:,3) = sq3*eta;

  if(method_order>=3) % 3rd order
    phi(:,:,4) =  3.0*xi.*eta;
    phi(:,:,5) =  sq5*(1.5*xi2 - 0.5);
    phi(:,:,6) =  sq5*(1.5*eta2 - 0.5);

  if(method_order>=4) % 4th order
    phi(:,:,7 )= sq3*sq5*eta.*(1.5*xi2 - 0.5);
    phi(:,:,8 )= sq3*sq5*xi.*(1.5*eta2 - 0.5);
    phi(:,:,9 )= sq7*(2.5*xi3 - 1.5*xi);
    phi(:,:,10)= sq7*(2.5*eta3 - 1.5*eta);

  if(method_order>=5) % 5th order
    phi(:,:,11)= sq3*sq7*(2.5*xi3 - 1.5*xi).*eta;
    phi(:,:,12)= sq3*sq7*(2.5*eta3 - 1.5*eta).*xi;
    phi(:,:,13)= 5.0/4.0*(3.0*xi2 - 1.0).*(3.0*eta2 - 1.0);
    phi(:,:,14)= 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0;
    phi(:,:,15)= 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0;

  end; % 5
  end; % 4
  end; % 3
  end; % 2

end
