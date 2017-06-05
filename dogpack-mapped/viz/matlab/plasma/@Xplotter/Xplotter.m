% I am using Xplotter as a dummy class
% to package the routines that are used to
% plot fluxes and the time history of the Xpoint.
%
% This is just a mechanism to allow two different user-accessible
% routines to access the same subroutines without exposing the
% subroutines to everyone.
%
function xplotter = Xplotter()
  xplotter = struct;
  xplotter = class(xplotter,'Xplotter');
end
