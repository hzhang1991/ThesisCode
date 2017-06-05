
function time_label = get_time_label(xplotter,do_rescale_times)
  time_label='time';
  if(do_rescale_times)
    % time_label='time*\Omega_i';
    %time_label='t*\Omega_i (time * ion gyrofrequency)';
    time_label='time per gyroperiod';
  end
end

