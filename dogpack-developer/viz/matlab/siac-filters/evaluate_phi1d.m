function pl = evaluate_phi1d(mpp, xg)
%EVALUATE_PHI.  Returns the basis functions evaluated at each quadrature point.
%
% This function evaluates the basis functions at each point specified by xg. It
% is assumed that xg \in [-1,1], the canonical interval.
%
% Input:
% ------
%
%     mpp - degree of polynomial.  mpp \in {0, 1, ..., whatever }.
%     xg  - Quadrature points to evaluate the polynomial
%
% Returns:
% --------
%
% pl( mpp+1, 6 ) - Legendre polynomials evaluated at each point in xg.
%
% See also: gleg_points.m

pl = zeros(mpp+2,6);
pl(1,:)=1.d0;
if (mpp >= 1) 
   pl(2,:)=xg;
end
if (mpp >= 2)
   for j=1:6
     for m=1:mpp
        pl(m+2,j)=((2*m+1.d0)*pl(m+1,j)*xg(j)-m*pl(m,j))/(m+1);
     end
   end
end

end
