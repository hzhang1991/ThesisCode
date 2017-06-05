function LegVals = GetCart2Legendre(kmax, s2d)
%GETCART2LEGENDRE    Sample the 2D Cartesian Legendre polynomials.
%
% LegVals = GETCART2LEGENDRE( kmax, s2d ) samples the 2D Legendre polynomials
% used in DoGPack at the list of points provided by s2d.
%
% Input:
%
%    kmax - number of basis elements in the current element. This should agree with
%    method_order.  See also: get_method_order.
%
%    s2d - points to be sampled (written in canonical variables).
%
% Output:
%
%    LegVals(1:kmax, 1:numpts ) - the kth Legendre polynomial evaluated at
%    point m.
%
% See also: GetUnst2Legendre, SAMPLE_BASIS_FUNCTIONS2.

  method_order = get_method_order(kmax,2);

  LegVals = zeros(kmax, size(s2d,1) );
  
  sq3 = sqrt(3);
  sq5 = sqrt(5);
  sq7 = sqrt(7);

  % Canonical variables:
  xi  = s2d(:,1);
  eta = s2d(:,2);

  % Quick lookup for quadrature points:
  xi2 = xi .*xi;
  xi3 = xi2.*xi;
  xi4 = xi3.*xi;

  eta2 =  eta.*eta;
  eta3 = eta2.*eta;
  eta4 = eta3.*eta;

  switch( method_order )

      case 1
          LegVals(1,:) = 1.0 + 0*xi;

      case 2
        LegVals(1,:)  = 1.0;      
        LegVals(2,:)  = sq3*xi;
        LegVals(3,:)  = sq3*eta;      

      case 3
        LegVals(1,:)  = 1.0;      
        LegVals(2,:)  = sq3*xi;
        LegVals(3,:)  = sq3*eta;      
        LegVals(4,:)  = 3.0*xi.*eta; 
        LegVals(5,:)  = sq5*(1.5*xi2 - 0.5);
        LegVals(6,:)  = sq5*(1.5*eta2 - 0.5);
      case 4

        LegVals(1,:)  = 1.0;      
        LegVals(2,:)  = sq3*xi;
        LegVals(3,:)  = sq3*eta;      
        LegVals(4,:)  = 3.0*xi.*eta; 
        LegVals(5,:)  = sq5*(1.5*xi2 - 0.5);
        LegVals(6,:)  = sq5*(1.5*eta2 - 0.5);
        LegVals(7,:)  = sq3*sq5*eta.*(1.5*xi2 - 0.5);      
        LegVals(8,:)  = sq3*sq5*xi.*(1.5*eta2 - 0.5);      
        LegVals(9,:)  = sq7*(2.5*xi3 - 1.5*xi);
        LegVals(10,:) = sq7*(2.5*eta3 - 1.5*eta);

      case 5

        LegVals(1,:)  = 1.0;      
        LegVals(2,:)  = sq3*xi;
        LegVals(3,:)  = sq3*eta;      
        LegVals(4,:)  = 3.0*xi.*eta; 
        LegVals(5,:)  = sq5*(1.5*xi2 - 0.5);
        LegVals(6,:)  = sq5*(1.5*eta2 - 0.5);
        LegVals(7,:)  = sq3*sq5*eta.*(1.5*xi2 - 0.5);      
        LegVals(8,:)  = sq3*sq5*xi.*(1.5*eta2 - 0.5);      
        LegVals(9,:)  = sq7*(2.5*xi3 - 1.5*xi);
        LegVals(10,:) = sq7*(2.5*eta3 - 1.5*eta);
        LegVals(11,:) = sq3*sq7*(2.5*xi3 - 1.5*xi).*eta;
        LegVals(12,:) = sq3*sq7*(2.5*eta3 - 1.5*eta).*xi;
        LegVals(13,:) = 5.0/4.0*(3.0*xi2 - 1.0).*(3.0*eta2 - 1.0);      
        LegVals(14,:) = 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0;      
        LegVals(15,:) = 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0;

  end

end % End of function
