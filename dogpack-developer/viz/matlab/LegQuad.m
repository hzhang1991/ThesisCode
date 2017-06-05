function [y,w] = LegQuad(N,a,b,tol,maxiter)
%LEGQUAD.   Compute the Legendre-Gauss quadrature points.
%
% Compute the Legendre-Gauss quadrature points of order N with
%         w(x) = 1 ,    x in (-1,1).
% The resulting quadrature points are then mapped to y in [a,b].
%
% With only one input value (N), the [a,b] default to the interval [-1,1].
%
% Usage: [y,w] = LegQuad(N,a,b,tol,maxiter)
%
% Inputs:
%
%    N       - number of qudrature points
%    [a,b]   - interval of integration.  Default [a,b] = [-1,1]
%    tol     - tolerance used for Newton iteration
%    maxiter - maximum number of iterations used in Newton iteration.
%    
% Output:
%
%    y - quadrature points
%    w - quadrature weights.
%
% Integration is then accomplished by \sum_i f( y_i ) w_i.
%
%----------------------Description of the algorithm.----------------------%
%
% The zeros of the Nth Legendre Polynomial P_{N}(x) are computed using:
%
%     P_{n+1} = [(2n+1)/(n+1)]x P_{n} - [n/(n+1)] P_{n-1}.
%
%
% Using a very accurate initial guess, the zeros are found by standard
% application of Newton's method:
%
%                   xnew = xold - P_{N}(xold) / P_{N}'(xold)
%
% The derivative P_{N}'(x) is also found using recursion:
%
%   P_{N}' = N/(1-x^2) ( -x P_{N} + P_{N-1} )
%
% Weights are then obtained by evaluating
%
%               w_i = c_i / ( (1-x^2)*(P_{N}'(x))^2 )
%
% at x = x_i.  Instead of using c_i, the weights are normalized according to
%
%     sum{w_i} =  \int_{-1}^{1} dx = 2.
%
% Author: Mathew F. Causley, Michigan State University Fall 2012.
%
%-------------------------------------------------------------------------%

% Set default values for unspecified parameters
if nargin<2,    a  = -1;    b  = 1;     end
if nargin<4,    tol= eps;               end
if nargin<5,    maxiter= 5;             end

%-------------------Initial guess is O(N^-4) accurate.--------------------%
%
% THIS IS THE PART THAT IS SPECIAL ABOUT THIS CODE.
%
%-------------------------------------------------------------------------%
th  = (4*N-1:-4:3)'*pi/(4*N+2);
x   = cos( th+2*cot(th)/(4*N+2)^2 );	%initial guess for the zeros/nodes

%------------Use recursion to evaluate P_{N}(x) and P_{N}'(x)-------------%
for k =1:maxiter
	Pnm1= ones(N,1);	%P_0 = 1
	Pn  = x;			%P_1 = x
	for n=1:N-1  %compute P_{n+1}(x) at the approximate roots x=x_i.
		Pnp1= ( (2*n+1)*x.*Pn- n*Pnm1 )/(n+1);
		Pnm1= Pn;
		Pn  = Pnp1;
	end
	%--------Construct Q = P_{N}' and update x with Newton's method-------%
	Q  = N*(-x.*Pn + Pnm1)./( (1-x.^2) );
	x  = x-Pn./Q;
    
	if norm(Pn./Q,inf)<tol
		break;
	end
end

%--------------Map from [-1,1] to [a,b], and compute weights--------------%
y = (1-x)*a/2+(1+x)*b/2;
w = 1./((1-x'.^2).*Q'.^2);
w = w/sum(w)*(b-a);	%Normalize weights

end
