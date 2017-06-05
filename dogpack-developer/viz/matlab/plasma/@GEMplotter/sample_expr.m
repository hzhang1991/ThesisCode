function [value,time,var_name,display_string] ...
  = sample_expr(c,eval_string,frame_number,varargin)

  [state,time] = read_state(c,frame_number);

  [value,var_name,display_string] ...
    = sample_expr_from_state(c,eval_string,state,varargin{:});
end
