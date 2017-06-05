
function [var_vals,time]  = sample(c, var_name, frame_number, varargin)
  [state,time] = read_state(c,frame_number);
  var_vals  = sample_var_from_state(c, var_name, state, varargin{:});
end
