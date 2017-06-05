function q0 = q_init1d( x )
%Q_INIT.  Initial conditions
%
% This function can be used to compute exact solutions. 
%

    q0 = sin( 2.0*pi*x );
    % q0 = exp( sin( 2.0*pi*x ) );

%   width = 2.0*0.25;
%   s     = mod(x-0.5,1)

%   I = find( abs(s) > width/2.0 );
%   q0 = 0.*x;
%   q0(I) = cos(pi*s(I)/width).^6;

end
