function LegVals=GetCart3Legendre(kmax, s3d)
%GETCART3LEGENDRE    Sample the 3D Legendre polynomials.
%
% LegVals = GETCART3LEGENDRE( kmax, s3d ) samples the 3D Cartesian 
% Legendre polynomials used in DoGPack at the list of points provided by s3d.
%
% Input:
%
%    kmax - number of basis elements in the current element. This should agree with
%    method_order.  See also: get_method_order.
%
%    s3d - points to be sampled (written in canonical variables).  
%          size(s3d) = [numpts, 3];
%
% Output:
%
%    LegVals(1:kmax, 1:numpts ) - the kth Legendre polynomial evaluated at
%    point m.


  meth1 = get_method_order(kmax,3);
  numpts  = size( s3d, 1 );
  LegVals = zeros(kmax, numpts );
  
  sq3 = sqrt(3);
  sq5 = sqrt(5);
  sq7 = sqrt(7);
  
  for m=1:(numpts)
    
    xi   = s3d(m,1);
    eta  = s3d(m,2);
    zeta = s3d(m,3);
    
    xi2 = xi*xi;
    xi3 = xi2*xi;
    xi4 = xi3*xi;
    
    eta2 = eta*eta;
    eta3 = eta2*eta;
    eta4 = eta3*eta;
    
    zeta2 = zeta*zeta;
    zeta3 = zeta2*zeta;
    zeta4 = zeta3*zeta;
    
    onehalf = 0.5;
    
    switch meth1
     case 4
      LegVals(1,m)  = 1.0;      
      LegVals(2,m)  = sq3*xi;
      LegVals(3,m)  = sq3*eta;      
      LegVals(4,m)  = sq3*zeta;
      LegVals(5,m)  = 3.0*xi*eta;
      LegVals(6,m)  = 3.0*xi*zeta;
      LegVals(7,m)  = 3.0*eta*zeta;
      LegVals(8,m)  = onehalf*sq5*(3.0*xi2-1.0);
      LegVals(9,m)  = onehalf*sq5*(3.0*eta2-1.0);
      LegVals(10,m) = onehalf*sq5*(3.0*zeta2-1.0);
      LegVals(11,m) = onehalf*sq3*sq5*eta*(3.0*xi2-1.0);
      LegVals(12,m) = onehalf*sq3*sq5*zeta*(3.0*xi2-1.0);
      LegVals(13,m) = onehalf*sq3*sq5*xi*(3.0*eta2-1.0);
      LegVals(14,m) = onehalf*sq3*sq5*zeta*(3.0*eta2-1.0);
      LegVals(15,m) = onehalf*sq3*sq5*xi*(3.0*zeta2-1.0);
      LegVals(16,m) = onehalf*sq3*sq5*eta*(3.0*zeta2-1.0);
      LegVals(17,m) = 3.0*sq3*xi*eta*zeta;
      LegVals(18,m) = onehalf*sq7*xi*(5.0*xi2-3.0);
      LegVals(19,m) = onehalf*sq7*eta*(5.0*eta2-3.0);
      LegVals(20,m) = onehalf*sq7*zeta*(5.0*zeta2-3.0);
     case 3
      LegVals(1,m)  = 1.0;      
      LegVals(2,m)  = sq3*xi;
      LegVals(3,m)  = sq3*eta;      
      LegVals(4,m)  = sq3*zeta;
      LegVals(5,m)  = 3.0*xi*eta;
      LegVals(6,m)  = 3.0*xi*zeta;
      LegVals(7,m)  = 3.0*eta*zeta;
      LegVals(8,m)  = onehalf*sq5*(3.0*xi2-1.0);
      LegVals(9,m)  = onehalf*sq5*(3.0*eta2-1.0);
      LegVals(10,m) = onehalf*sq5*(3.0*zeta2-1.0);
     case 2
      LegVals(1,m)  = 1.0;      
      LegVals(2,m)  = sq3*xi;
      LegVals(3,m)  = sq3*eta;      
      LegVals(4,m)  = sq3*zeta;
     case 1
      LegVals(1,m) = 1.0;
    end
  end

end
