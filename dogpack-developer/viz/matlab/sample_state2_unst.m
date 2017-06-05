function qvals = sample_state2_unst(q, method_order, kmax, LegVals)
%SAMPLE_STATE_UNST   Sample the basis functions.
%
%QVALS = SAMPLE_STATE2_UNST(q, method_order, kmax, LegVals )
%
% Turn Galerkin coefficients into point values on an unstructured grid
%
% Parameters:
%
%   q( melems*mpts, 1:meqn(or maux), 1:kmax ) - the basis functions
%   method_order - order of the method (this should be consistent with kmax)
%   kmax         - number of basis functions per element
%   LegVals      - Legendre polynomials evaluated at a list of points
%
% Returns:
%
%   qvals( melems*mpts, meqn ) - the function evaluated at all melems*mpts
%                                points, for each equation meqn (or maux).
%

   melems = size(q,1);
   meqn   = size(q,2);
   mpts   = size(LegVals,2);
   
   points_per_dir = sqrt(mpts);
   
   qvals = zeros(melems*mpts, meqn);

%TODO - it might be nice to vectorize this process.  It should speed up the
%       code, but we need to perform multi-d matrix multiplication to do so.
%       (-DS)

%  index = 0;
%  for i=1:melems
%    for m=1:mpts
%      index = index+1;
%      for j=1:meqn
%        
%        qvals(index,j) = 0;
%        for k=1:kmax
%          qvals(index,j) = qvals(index,j) + q(i,j,k)*LegVals(k,m);
%        end  
%        
%      end
%    end
%  end

  % Vectorized version of the above section that has three for loops.  I'm
  % assuming that the for loop over meqn is small compared to the other
  % directions.
  %
  % TODO: THIS SHOULD BE CHECKED ON A SYSTEM WITH MORE THAN ONE EQUATION.
  % (-DS)
  qm = zeros( melems, kmax );
  for me=1:meqn

    % Sample the 'mth' equation:
    qm(:,:)  = squeeze( q(:,me,:) );
    qm_pts   = qm * LegVals;
    qm_pts   = reshape( qm_pts', melems*mpts, 1 );

    % Save the points from the 'mth' equation:
    qvals( :, me) = qm_pts(:);

  end

end
