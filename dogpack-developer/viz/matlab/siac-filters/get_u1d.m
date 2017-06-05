function [u, time] = get_u1d( fname, Elmts, morder )
%GET_U.  Function to read u from a file.
%
% This function was written to be able to interface with DoGPack, where the
% polynomials are rescaled to be orthonormal.  Additionally, DoGPack saves the
% final time of the solution.
%
% Parameters
% ----------
%
%           fname - name of file to extract data from
%
%           Elmts - number of grid elements for that solution.  (This is used
%           to reshape the vector of unkowns.
%
%           morder - order of accuracy.  (morder == mpp+1)
%
% Returns
% -------
%
%           u - vector of unkowns. u = u( 1:morder, 1:Elmts, meqn ).
%               NOTE: this is a vector of (nonorthonormal) Legendre
%               coefficients!  They have been converted from what DoGPack
%               first constructs.
%
%       time  - time of solution.  DoGPack saves the time for each output, so
%               we may as well return that as well.
%
% See also: main_driver.

    fids = fopen( fname, 'r' );
    if fids==-1
        error(['File  ',fname,'  not found.']);
    end

    % Read in all of the appropriate data
    time = fscanf(fids,'%e',1);
    qtmp = fscanf(fids,'%e',[1,inf]);
    fclose(fids);

    qtmp = transpose(qtmp);
    % u    = reshape(qtmp, Elmts, morder, meqn);
    u    = reshape(qtmp, Elmts, morder);
    clear qtmp;

    % Rescale each polynomial in u
%   for me=1:meqn
    for k=1:morder
    for i=1:Elmts
            u(i,k) = u(i,k) * sqrt( 2.0*k-1.0 );
    end
    end
%   end

    % Add in ghost cells (periodic boundary conditions)
    u = [ u(Elmts,:); u; u(1,:) ];
    u = u';

end
