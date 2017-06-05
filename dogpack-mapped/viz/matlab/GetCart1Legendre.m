function phi = GetCart1Legendre(kmax, xi )
%GETCART1LEGENDRE    Sample 1D Legendre polynomials.
%
% phi = GETCART1LEGENDRE(kmax, xi ) samples the basis functions at a list of
% points provided in xi.
%
% Input:
%
%      xi = points that will be sampled
%    kmax = number of polynomials considered (In 1D, same as the order of the
%   method)
%
% Output:
%
%    phi(1:numpts, 1:kmax ) = array of the basis functions sampled at each
%    point.
%


  sq3 = sqrt(3);
  sq5 = sqrt(5);
  sq7 = sqrt(7);
  sq11 = sqrt(11);

  numpts = length(xi);
  phi    = zeros(numpts, kmax);

  xi2  = xi.*xi;
  xi3  = xi.*xi2;
  xi4  = xi.*xi3;
    
  phi(:,1) = 1;
    
  if( kmax>1 )
    phi(:,2)  = sq3*xi;
  end
    
  if (kmax>2)
    phi(:,3)  = sq5*(1.5*xi2 - 0.5);
  end
    
  if (kmax>3)
      phi(:,4)  = sq7*(2.5*xi3 - 1.5*xi);
  end
    
  if (kmax>4)
      phi(:,5) = 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0;
  end
    
  if (kmax>5)
      phi(:,6) = (63/8)*sq11 * ( xi.^5 - (10/9)*xi3 + (5/21)*xi );
  end

end
