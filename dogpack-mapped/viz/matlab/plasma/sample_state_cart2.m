function qvals = sample_state_cart2(q, how_sample, method_order)
%SAMPLE_STATE_CART2.  Function for sampling the Galerkin expansion.
%
% This function is an alternative to sample_state2_cart, but is only called
% from a single routine located in GEMCHALLENGE. These two function calls
% should be merged!  (-DS)
%
% qvals = SAMPLE_STATE_CART2(q, how_sample, method_order ) samples
% the Galerkin expansion at a colle
%
% how_sample:
% 0 == gaussian quadrature
% positive integer: points per cell
% negative integer: cells per point (not implemented)
%
% Input:
%
%    q( mx, my, kmax,  ) - what's the size of q when it gets passed in here?
%    Look at: read_state2_cart.m
%
%    how_sample - flag deciding gaussian or equispaced quadrature points.
%                 how_sample == 0 - gaussian points.
%                 how_sample ~= 0 - equispaced points.
%
%    method_order - order of the method.  (commonly called meth1 in the code).
%

    if(~exist('method_order','var'))

        % Are you sure you want to be using size(q) here?  I'm not planning on fixing
        % this if it's broken.  On top of it, I'm not sure what size q is
        % supposed to be when it gets passed in here ... (-DS)

        % function method_order = get_method_order(kmax, ndims)
        method_order = get_method_order( size(q,3) );
    else
        assert(method_order == get_method_order(size(q,3)));
    end
    if(~exist('how_sample','var'))
        how_sample = method_order;
    end

    if(how_sample == 0)
        % use gaussian quadrature points
        qvals = sample_state2(q, method_order, 2, method_order);
    else
        % use equispaced points
        sample_rate = how_sample;
        qvals = sample_state2(q, method_order, 1, sample_rate);
    end

end

function qvals = sample_state2(q, method_order, point_type, numpts)
%SAMPLE_STATE2.     Sample something ...
%
% point_type: 1=equispaced, 2=gaussian_quadrature
% numpts: if positive, the number of points per cell;
%         if negative, the number of cells per point;

  phi = sample_basis_functions2(method_order, point_type, numpts);
  mx=size(q,1);
  my=size(q,2);

  % create qvals from q
  %
  mx_out = numpts*mx;
  my_out = numpts*my;
  qvals=zeros(mx_out,my_out,size(q,4));
  % contract tensor product of phi with q.
  % perhaps I could use repmat here to avoid for loops,
  % but probably at the cost of greater memory usage,
  % and these are not large loops.
  for m1=1:numpts
  for m2=1:numpts
    qval=zeros(mx,my,1,size(q,4));
    %size(qval)
    %size(q)
    for k=1:get_kmax(method_order)
      qval=qval+phi(m1,m2,k)*q(:,:,k,:);
    end
    qvals(m1:numpts:end,m2:numpts:end,:) ...
      =reshape(qval,[mx,my,size(q,4)]);
  end
  end
end

function qfull = sample_state2_coarse(q, method_order, point_type, numpts)
% if numpts is odd this should sample at the center of the sampled cells
% but if numpts is even then this needs to sample at corners, which
% it seems best to do by averaging corner samples from four cells
% but which we do by just choosing one of the cells.

  error(['sample_state2_coarse not yet implemented']);
end

% I'm not clear that anyone ever called this function ... (-DS)
%
% We have renamed sample_state2_cart_mod to be sample_state2_cart in order to
% be consistent with the other 1 and 3d plotting routines.
%
%%%%function qvals = sample_state2_cart(q, method_order, point_type, numpts);
%%%%%
%%%%% Turn coefficients into point values on Cartesian grid
%%%%%
%%%%% point_type: 1=equispaced, 2=gaussian_quadrature
%%%%% numpts: if positive, the number of points per cell;
%%%%%         if negative, the number of cells per point;
%%%%%
%%%%% This function was moved from sample_state2_cart.m here, because the only
%%%%% application that was using it was located in the plasma visualization
%%%%% library. In order to clean up viz/matlab, I've moved it here.  How all of
%%%%% the sample_state[n]* routines match each other. (-DS)
%%%%%

%%%%  if(~exist('point_type'))
%%%%    point_type=1;
%%%%  end
%%%%  if(~exist('numpts'))
%%%%    numpts=method_order;
%%%%  end
%%%%  if(point_type==2)
%%%%    how_sample=0;
%%%%  else
%%%%    how_sample=numpts;
%%%%  end

%%%%  % call this directly instead
%%%%  qvals = sample_state_cart2(q, how_sample, method_order);

%%%%end
