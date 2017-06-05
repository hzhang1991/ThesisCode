function [wg, xg] = gleg_points()
%GL_POINTS.  Returns the Gauss Legendgre weights and points
%
% This routine is hard-coded to only construct a total of 6 quadrature points.
% The reason for this is we do not currently have the code to construct the
% evaluation of the SIAC kernel at arbitrary points, because those too have
% been hard coded in the past.
%
% Returns:
% --------
%
% wg( 6,1 ) - gauss Legendgre weights
% xg( 6,1 ) - gauss Legendgre points
%
% See also: evaluate_phi1d.

    wg=zeros(6,1);
    xg=zeros(6,1);

    xg(1)=-0.466234757101576013906150777246997304567d0*2;
    xg(2)=-0.330604693233132256830699797509952673503d0*2;
    xg(3)=-0.119309593041598454315250860840355967709d0*2;
    xg(4)=-xg(3);
    xg(5)=-xg(2);
    xg(6)=-xg(1);

    wg(1)=1.7132449237916972e-01/2.d0;
    wg(2)=3.6076157304813861e-01/2.d0;
    wg(3)=4.6791393457269098e-01/2.d0;
    wg(4)=wg(3);
    wg(5)=wg(2);
    wg(6)=wg(1);

end
