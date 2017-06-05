function q0 = q_init2d( x, y )
%Q_INIT.  Initial conditions
%
% This function can be used to compute exact solutions. 
%

    q0 = sin( 2.0*pi*x ) .* cos( 2.0*pi*y );

end
