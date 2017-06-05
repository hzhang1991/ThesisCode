
% plot a state variable
%
% example of how to use:
%
%  for i=1:20; plot(g,'rho_e',i); ui=input('prompt> '); if(ui==0) break; end; end
%
function plot(c,eval_string,frame_number,varargin)
  pause_time=-1;
  for k=1:numel(frame_number)
    % wait for user to hit return after showing first frame
    if(k>1)
      if(pause_time>=0) pause(pause_time);
      else
        prompt = ['view frame ' num2str(frame_number(k)) '? '];
        user_input = input(prompt, 's');
        if(numel(user_input)==0) user_input='y'; end
        first_char = user_input(1);
        if(first_char=='q' || first_char=='Q' || ...
           first_char=='n' || first_char=='N'); return;
        elseif(first_char=='y' || first_char=='Y');
        else
          [pause_time,status]=str2num(user_input);
          % do not flicker while showing movie
          %set(gcf,'DoubleBuffer','on');
          if(~status) error('invalid input'); end;
        end
      end
    end
    % [state,time] = read_state(c,frame_number(k));
    if(isa(eval_string,'cell'))
      c.s.stateInfo = initStateInfo(c.s.stateInfo);
      for display_idx=1:numel(eval_string)
        c.s.stateInfo = plot_helper(c,eval_string{display_idx},...
          frame_number(k), 'display_idx',display_idx,varargin{:});
      end
    else
      plot_helper(c,eval_string,frame_number(k),varargin{:});
    end
  end
end

% stateInfo is a mechanism to avoid repeatedly
% reading in the same data
%
function stateInfo = initStateInfo(stateInfo)
  stateInfo.state = [];
  stateInfo.time = -1;
  stateInfo.frame_number = -1;
  stateInfo.auxState = [];
end

function stateInfo = plot_helper(c,eval_string,frame_number,varargin)
  %profile on;
  s = c.s;
  % % no reason to do this since there is no speedup
  %varIdx_s = c.s.varMap.varIdx_s;
  %if(isfield(varIdx_s,eval_string))
  %  var_name=eval_string;
  %  [value,time] = sample(c, var_name, frame_number);
  %  display_string = c.s.plot_params.display_names{varIdx_s.(var_name)};
  %else
  % get sampled values
  c.s.stateInfo.frame_number = frame_number;
  [value,var_name,display_string,c.s.stateInfo] ...
    = sample_expr_from_state(c,eval_string,varargin{:});
  display_frame(c, c.s.stateInfo.time, value, var_name, 'display_string', ...
    display_string, varargin{:});
  stateInfo = c.s.stateInfo;
  %profile off;
  %profile viewer;
end

