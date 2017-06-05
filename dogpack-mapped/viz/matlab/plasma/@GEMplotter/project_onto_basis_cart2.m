
% quad_vals contains sample values at gaussian quadrature points
%
% state contains values of coefficients of the expansion in
%   legendre polynomials
%
function state = project_onto_legendre_basis(quad_vals,space_order)
  wght_2d = get_wght_2d(space_order);
  mx=size(quad_vals,1)/space_order;
  my=size(quad_vals,2)/space_order;
  kmax=get_kmax(space_order);
  meqns=size(quad_vals,3);
  phi = sample_basis_functions2(space_order, 2);
  for k_idx=1:kmax
    accumulator = zeros(mx,my,meqns);
    for n1=1:space_order
    for n2=1:space_order
      % The legendre functions are orthonormal.
      % So to get the coefficients of the legendre function expansion
      % we just integrate the state values multiplied by the legendre 
      % function values.  (We integrate by sampling at Gaussian
      % quadrature points and taking a linear combination with
      % Gaussian weights.)
      accumulator = accumulator + wght_2d(n1,n2)*phi(n1,n2,k_idx) ...
        *quad_vals(n1:space_order:end,n2:space_order:end,:);
    end
    end
    state(:,:,k_idx,:) = reshape(accumulator,[mx,my,1,meqns])/4.0;
  end
end

function wght_2d = get_wght_2d(space_order)
  wght_1d = get_wght_1d(space_order);
  wght_2d = transpose(wght_1d) * wght_1d;
end

function wght_1d = get_wght_1d(space_order)
    if(space_order==1)
      wght_1d = 2.0;
    elseif(space_order==2)
      wght_1d = [1., 1.];
    elseif(space_order==3)
      wght_1d = [5., 8., 5.]/9.0;
    elseif(space_order==4)
      w1 = (18.0 - sqrt(3)*sqrt(10))/36.0;
      w2 = (18.0 + sqrt(3)*sqrt(10))/36.0;
      wght_1d = [w1, w2, w2, w1]; 
    elseif(space_order==5)
      w1 = (322.0 - 13.0*sqrt(7)*sqrt(10))/900.0;
      w2 = (322.0 + 13.0*sqrt(7)*sqrt(10))/900.0;
      w3 = 128.0/225.0;
      wght_1d = [w1, w2, w3, w2, w1]; 
    else
      error(['invalid space_order: ' num2str(space_order)]);
    end
end

