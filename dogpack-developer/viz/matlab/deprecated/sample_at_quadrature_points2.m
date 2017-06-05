
function vals=sample_at_quadrature_points2(state)

  method_order = get_method_order(size(state,3))
  quad_vals=zeros(size(state,1),size(state,2),method_order,size(state,4));
  % 2d quadrature weights and points: tensor product
  %
  %weight = transpose(weight_1d) * weight_1d;
  %position = transpose(position_1d) * position_1d;
  %
  % evaluate Legendre basis function at each 2d quadrature point
  phi = sample_basis_functions2(method_order,2);
  %
  % contract tensor product of phi with q.
  for m1=1:method_order
  for m2=1:method_order
    qval=zeros(mx,my,1,size(q,4));
    for k=1:get_kmax(method_order)
      qval=qval+phi(m1,m2,k)*q(:,:,k,:);
    end
    %qfull(m1:method_order:end,m2:method_order:end,:) ...
    %  =reshape(qval,[mx,my,size(q,4)]);
  end
  end
  %
end

% should reimplement sample_basis_functions in sample_state2 using
% this function to avoid code redundancy
%
function phi = sample_basis_functions_at_points(position_1d, method_order)
  %
  assert(size(position_1d)==method_order);
  %
  %
end
