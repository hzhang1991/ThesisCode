function LegVals = GetUnst2Legendre(kmax, z)
%GETUNST2LEGENDRE   Sample the 2D unstructured Legendre polynomials
%
%LEGVALS = GETUNST2LEGENDRE(kmax, z )
%
%  Sample 2D Legendre polynomials that are on an unstructured grid.
%
% Parameters:
%
%    kmax                      - number of polynomails sampled
%    z(1:numpts, 2)            - points that are being sampled.
%                     
%
% Returns:
%
%    LegVals(1:kmax, 1:numpts) - The legendre polynomials evaluated 
%                                at every point indicated in the set of 
%                                points listed by z.
%
% See also: GetCart2Legendre
 
  meth1 = get_method_order(kmax,2);
  numpts  = size( z, 1 );
  MonVals = zeros(kmax, numpts);

% for m=1:(points_per_dir^2)
%   
%   xi  = z(m,1);
%   eta = z(m,2);
%   
%   switch meth1
%    case 5
%     MonVals(1,m)  = 1.0;
%     MonVals(2,m)  = xi;
%     MonVals(3,m)  = eta;
%     MonVals(4,m)  = xi*eta;
%     MonVals(5,m)  = xi^2;
%     MonVals(6,m)  = eta^2;
%     MonVals(7,m)  = eta*xi^2;
%     MonVals(8,m)  = xi*eta^2;
%     MonVals(9,m)  = xi^3;
%     MonVals(10,m) = eta^3;
%     MonVals(11,m) = xi^3*eta;
%     MonVals(12,m) = eta^3*xi;
%     MonVals(13,m) = xi^2*eta^2;
%     MonVals(14,m) = xi^4;
%     MonVals(15,m) = eta^4;
%    case 4
%     MonVals(1,m)  = 1.0;
%     MonVals(2,m)  = xi;
%     MonVals(3,m)  = eta;
%     MonVals(4,m)  = xi*eta;
%     MonVals(5,m)  = xi^2;
%     MonVals(6,m)  = eta^2;
%     MonVals(7,m)  = eta*xi^2;
%     MonVals(8,m)  = xi*eta^2;
%     MonVals(9,m)  = xi^3;
%     MonVals(10,m) = eta^3;
%    case 3
%     MonVals(1,m)  = 1.0;
%     MonVals(2,m)  = xi;
%     MonVals(3,m)  = eta;
%     MonVals(4,m)  = xi*eta;
%     MonVals(5,m)  = xi^2;
%     MonVals(6,m)  = eta^2;      
%    case 2
%     MonVals(1,m)  = 1.0;
%     MonVals(2,m)  = xi;
%     MonVals(3,m)  = eta;
%    case 1
%     MonVals(1,m) = 1.0;
%   end
% end

  % Vectorized version of what's written above:
  xi  = z(:, 1);
  eta = z(:, 2);
  switch meth1

     case 5
      MonVals(1,:)  = 1.0;
      MonVals(2,:)  = xi(:);
      MonVals(3,:)  = eta(:);
      MonVals(4,:)  = xi(:).*eta(:);
      MonVals(5,:)  = xi(:).^2;
      MonVals(6,:)  = eta(:).^2;
      MonVals(7,:)  = eta(:).*xi(:).^2;
      MonVals(8,:)  = xi(:).*eta(:).^2;
      MonVals(9,:)  = xi(:).^3;
      MonVals(10,:) = eta(:).^3;
      MonVals(11,:) = xi(:).^3.*eta(:);
      MonVals(12,:) = eta(:).^3.*xi(:);
      MonVals(13,:) = xi(:).^2.*eta(:).^2;
      MonVals(14,:) = xi(:).^4;
      MonVals(15,:) = eta(:).^4;
     case 4
      MonVals(1,:)  = 1.0;
      MonVals(2,:)  = xi(:);
      MonVals(3,:)  = eta(:);
      MonVals(4,:)  = xi(:).*eta(:);
      MonVals(5,:)  = xi(:).^2;
      MonVals(6,:)  = eta(:).^2;
      MonVals(7,:)  = eta(:).*xi(:).^2;
      MonVals(8,:)  = xi(:).*eta(:).^2;
      MonVals(9,:)  = xi(:).^3;
      MonVals(10,:) = eta(:).^3;
     case 3
      MonVals(1,:)  = 1.0;
      MonVals(2,:)  = xi(:);
      MonVals(3,:)  = eta(:);
      MonVals(4,:)  = xi(:).*eta(:);
      MonVals(5,:)  = xi(:).^2;
      MonVals(6,:)  = eta(:).^2;      
     case 2
      MonVals(1,:)  = 1.0;
      MonVals(2,:)  = xi(:);
      MonVals(3,:)  = eta(:);
     case 1
      MonVals(1,:) = 1.0;
    end

  % Get transformation matrix to map monomial
  % values to Legendre values
  Mmat = GetMonomialToLegendre();
  Mmat = Mmat(1:kmax,1:kmax);

  % Transform the monomial values to legendre values:
  LegVals = zeros(kmax, numpts );
  for m=1:(numpts)
    LegVals(:,m) = Mmat*MonVals(:,m);
  end  

end

function Mmat = GetMonomialToLegendre()
%GETMONOMIALTOLEGENDRE.   Get Monomial to Legendre conversion matrices.
%
% Calling format: Mmat = GetMonomialToLegendre()
%
% This function is used for the 2D unstructured plotting routines.  It returns
% the transformation map necessary to convert monomials to legendre basis
% functions.
%
% The monomial basis functions (1, xi, eta, xi*eta, xi^2, eta^2, ... ), are 
% not orthogonal over a single element, but are easier to list and sample.
%
% The mapping Mmat : \R^{kmax x kmax} \to \R^{ kmax x kmax }
%
% converts the monomials to legendre basis elements.
%
% See also: GetUnst2Legendre, get_kmax

  % Quick lookup tables:
  sq2 = sqrt(2.0);
  sq3 = sqrt(3.0);
  sq5 = sqrt(5.0);
  sq7 = sqrt(7.0);
  sq10 = sqrt(10.0);
  sq13 = sqrt(13.0);
  sq19 = sqrt(19.0);
  sq23 = sqrt(23.0);
  sq71 = sqrt(71.0);

  % The following mapping was precomputed.
  % This can be truncated later by looking at Mmat(1:kmax, 1:kmax).
  % Here, we assume that kmax <= get_kmax( space_order, 2 ).
  Mmat = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          0.0, 3.0*sq2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          0.0, sq2*sq3, 2.0*sq2*sq3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          (5.0/21.0)*sq7, 8.0/sq7, 8.0/sq7, 60.0/sq7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          -10.0/9.0*(sq3/sq7), -2.0/(sq3*sq7), 4.0*sq3/sq7, 30.0*sq3/sq7, 5.0*sq3*sq7, 0.0, 0.0, 0.0, ...
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          -2.0/9.0*sq3*sq5, 2.0*sq5/sq3, 0, 6.0*sq3*sq5, sq3*sq5, ...
          6.0*sq3*sq5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          -8.0*sq3*sq13/117.0, 20.0/(sq3*sq13), -20.0/(sq3*sq13), 20.0*sq3/sq13, 40.0*sq3/sq13, 0.0, ...
          210.0*sq3/sq13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          -8.0/351.0*sq3*sq5*sq13*sq19, -4.0/247.0*sq3*sq5*sq13*sq19, 4.0/247.0*sq3*sq5*sq13*sq19, ...
          (20.0/39.0)*sq3*sq5*sq13*sq19, (344.0/741.0)*sq3*sq5*sq13*sq19, ...
          (32.0/57.0)*sq3*sq5*sq13*sq19, (602.0/247.0)*sq3*sq5*sq13*sq19, (56.0/19.0)*sq3*sq5*sq13*sq19, ...
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          0.0, -52.0/(sq3*sq19), -(8.0/19.0)*sq3*sq19, (8.0/3.0)*sq3*sq19, -4.0/(sq3*sq19), ...
          80.0/(sq3*sq19), (392.0/19.0)*sq3*sq19, (140.0/19.0)*sq3*sq19, 14.0*sq3*sq19, 0.0, ...
          0.0, 0.0, 0.0, 0.0, 0.0;
          %
          0.0, -(4.0/3.0)*sq7, -(8.0/3.0)*sq7, 8.0*sq7, 4.0*sq7, 0, 24.0*sq7, 60.0*sq7, 2.0*sq7, ...
          40.0*sq7, 0.0, 0.0, 0.0, 0.0, 0.0;
          %
          -(2.0/27.0)*sq5*sq7, -(26.0/21.0)*sq5*sq7, -(4.0/21.0)*sq5*sq7, -6.0*sq5*sq7, sq5*sq7, ...
          0.0, 0.0, 0.0, 16.0*sq5*sq7, 0, 72.0*sq5*sq7, 0.0, 0.0, 0.0, 0.0;
          %
          -(2.0/1053.0)*sq5*sq7*sq13*sq23, (314.0/6279.0)*sq5*sq7*sq13*sq23, ...
          -(544.0/6279.0)*sq5*sq7*sq13*sq23, -(2.0/13.0)*sq5*sq7*sq13*sq23, ...
          -(47.0/897.0)*sq5*sq7*sq13*sq23, (70.0/897.0)*sq5*sq7*sq13*sq23, 0.0, 0.0, ...
          -(752.0/897.0)*sq5*sq7*sq13*sq23, (1120.0/897.0)*sq5*sq7*sq13*sq23, ...
          -(1128.0/299.0)*sq5*sq7*sq13*sq23, (1680.0/299.0)*sq5*sq7*sq13*sq23, 0.0, 0.0, 0.0;
          %
          (55.0/74763.0)*sq7*sq13*sq71, -(440.0/2769.0)*sq7*sq13*sq71, -(440.0/2769.0)*sq7*sq13*sq71, ...
          -(140.0/923.0)*sq7*sq13*sq71, -(160.0/2769.0)*sq7*sq13*sq71, -(160.0/2769.0)*sq7*sq13*sq71, ...
          (360.0/71.0)*sq7*sq13*sq71, (360.0/71.0)*sq7*sq13*sq71, ...
          (6800.0/2769.0)*sq7*sq13*sq71, (6800.0/2769.0)*sq7*sq13*sq71, ...
          (10200.0/923.0)*sq7*sq13*sq71, (10200.0/923.0)*sq7*sq13*sq71, (1620.0/71.0)*sq7*sq13*sq71, 0.0, 0.0;
          %
          (2.0/639.0)*sq3*sq5*sq23*sq71, (4.0/14697.0)*sq3*sq5*sq23*sq71, ...
          -(280.0/14697.0)*sq3*sq5*sq23*sq71, -(20.0/71.0)*sq3*sq5*sq23*sq71, ...
          -(466.0/1633.0)*sq3*sq5*sq23*sq71, -(40.0/1633.0)*sq3*sq5*sq23*sq71, ...
          (60.0/71.0)*sq3*sq5*sq23*sq71, (60.0/sq71)*sq3*sq5*sq23, ...
          -(4.0/1633.0)*sq3*sq5*sq23*sq71, (280.0/1633.0)*sq3*sq5*sq23*sq71, ...
          (9780.0/1633.0)*sq3*sq5*sq23*sq71, (1260.0/1633.0)*sq3*sq5*sq23*sq71, ...
          (270.0/sq71)*sq3*sq5*sq23, 3.0*sq3*sq5*sq23*sq71, 0.0;
          %
          (2.0/9.0)*sq5, -(4.0/3.0)*sq5, 0.0, -20.0*sq5, -2.0*sq5, -20.0*sq5, 60.0*sq5, ...
          60.0*sq5, 12.0*sq5, 0.0, 60.0*sq5, 420.0*sq5, 270.0*sq5, 3.0*sq5, 210.0*sq5];

end
