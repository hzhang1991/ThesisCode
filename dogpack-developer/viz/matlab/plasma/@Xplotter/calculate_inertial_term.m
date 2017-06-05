
function mu3i_t = calculate_inertial_term(xplotter, t, mu3i);
  mu3i_t_middle = (mu3i(3:end)-mu3i(1:end-2)) ...
              ./ (  t(3:end)-  t(1:end-2));
  % extrapolate for boundary values
  mu3i_t_left = 2*mu3i_t_middle(1) - mu3i_t_middle(2);
  mu3i_t_rght = 2*mu3i_t_middle(end) - mu3i_t_middle(end-1);
  mu3i_t = zeros(size(mu3i));
  mu3i_t(1) = mu3i_t_left;
  mu3i_t(end) = mu3i_t_rght;
  mu3i_t(2:end-1) = mu3i_t_middle;
end
  
