function GenerateLegQuad( ptslist, tol )
%GENERATELEGQUAD.    Generate Legendre quadrature rules
%
% This routine can be used to generate arbitrary order Quadrature rules.
% The output is designed to generate code that can be integrated with DoGPack.
%
% See $DOGPACK/lib/Quadrature.cpp.
%
% Input:
% ------
%
%   ptslist : a list of points used to generate quadrature rules.
%   tol     : tolerance used for numerically finding the quadrature rules.
%             Default: 1e-14.
%
% Returns:
% --------
%
%   Nothing, save a bunch of stuff printed to the screen.
%
% See also: LegQuad.

    if( nargin > 1 )
        tol = 1e-14;
    end

    for m = 1:length( ptslist )

        numpts = ptslist(m);

        [y,w] = LegQuad(numpts, -1, 1, 1e-14, 10000 );

        fprintf(1, [...
'        case ', num2str( numpts, '%d' ), ':\n' ] );

        for n=1:numpts
            fprintf(1, [...
'            x1d.set(', num2str(n,'%d'), ', ', num2str(y(n), '%2.15e' ), ');\n'] );
        end
        fprintf(1,'\n');
        for n=1:numpts
            fprintf(1, [...
'            s1d.set(', num2str(n,'%d'), ', ', num2str(w(n), '%2.15e' ), ');\n'] );
        end

        fprintf(1, [...
'            break;\n\n'] );

end
