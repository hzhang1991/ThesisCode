function h=pdecont_tmp(p,t,u,lev)
%PDECONT Contour plot of PDE solution.
%
%       PDECONT(P,T,U) draws 10 level curves of the function represented
%       by the column vector U. P and T are point and triangles of a
%       mesh. The function plotted is extended from the values at the
%       points by linear interpolation.
%
%       PDECONT(P,T,U,N) plots using N levels.
%
%       PDECONT(P,T,U,V) plots using the levels specified by V.
%
%       H=PDECONT(P,T,U) returns handles to the line objects.





if nargin<4
  lev=10;
end

hh=pdeplot(p,[],t,'xydata',u,'xystyle','off','contour','on', ...
           'levels',lev,'colorbar','off');

if nargout>0
  h=hh;
end

