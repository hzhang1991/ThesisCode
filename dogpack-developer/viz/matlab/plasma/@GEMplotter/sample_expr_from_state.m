% if eval_string is simply the name of a variable
% we could handle this simple case separately,
% but there is no real reason to do so since we just
% end up evaluating an assignment, which incurs negligible
% cost since matlab does copy-on-write.
% 
% return everything needed for display.
function [value,var_name,display_string,stateInfo] = sample_expr_from_state(...
  c,eval_string,varargin)

  stateInfo = c.s.stateInfo;
  % parse keyword arguments
  %
  idx=1;
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      case 'var_name'; idx=idx+1;
        var_name=varargin{idx};
      otherwise
        %error(['invalid keyword: ' varargin{idx}]);
        % unrecognized keyword arguments are presumed to have one argument
        % and are ignored
        idx=idx+1;
    end
    idx=idx+1;
  end

  %var_vals = cell(size(var_names));

  % extract variables names from user-defined expression
  % and execute eval statements to ensure that the variables
  % referred to are defined in the context in which the
  % user-defined expression will be evaluated.
  %
  % we should also extract generic names like aux2 and q3
  % and use generic code to define such variables
  %
  [var_names, first_name] = extract_var_names(c, eval_string);
  if(~exist('var_name','var')); var_name=first_name; end
  params = c.s.params;
  m_e = params.elc_mass;
  m_i = params.ion_mass;
  % if numel(aux_names)>0
  %   astate = read_aux_state(c, frame_number);
  %   for i=1:numel(aux_names)
  %     aux_num = get_aux_num(aux_names{i});
  %     command_str=[aux_names{i} '= sample_aux(c, astate, aux_num);'];
  %     eval(command_str);
  %   end
  % end
  for i=1:numel(var_names)
    if(isfield(c.s.auxVarMap.varIdx_s,var_names{i}))
      command_str=['[' var_names{i} ...
        ', stateInfo] = sample_aux(c, var_names{i});'];
    else
      if(isempty(stateInfo.state))
        [stateInfo.state,stateInfo.time] = read_state(c,stateInfo.frame_number);
      end
      command_str=[var_names{i} ...
        ' = sample_var_from_state(c, var_names{i}, stateInfo.state);'];
    end
    eval(command_str);
    %[var_vals{i}, time] = sample(c, var_names{i}, frame_number);
  end
  %
  % evaluate user-supplied expression
  %
  value = eval(eval_string);
  %
  % replace variable names with descriptions in display string
  %
  display_string = replace_var_names(c,var_names,eval_string);
end

function [var_vals, stateInfo] = sample_aux(c, var_name)
    stateInfo = c.s.stateInfo;
    frame_number = stateInfo.frame_number;
    if(isempty(stateInfo.auxState))
      [stateInfo.auxState, stateInfo.time] = read_aux_state(c, frame_number);
    end
    components = c.s.auxVarMap.stateIndices_s.(var_name);
    var_coef = stateInfo.auxState(:,:,:,components);
    var_vals = sample_state_cart2(var_coef, c.s.sample_rate);
end

% function [component_vals, stateInfo] = sample_aux(c, varname);
%   component_state = astate(:,:,:,auxnum);
%   component_vals = sample_state_cart2(component_state, c.s.sample_rate);
% end

function all_names = extract_names(c, eval_string)
  name_arr = regexp(eval_string,'([A-Za-z]\w*)','tokens');
  name_arr = [name_arr{:}]; % convert cells to strings
  all_names = unique(name_arr); % sort -u
end

function [astate,time] = read_aux_state(c, frame_number)
  params = c.s.params;
  components = 1:params.maux;
  [astate,time]=read_state2_cart(params.datafmt, params.outputdir, ...
    frame_number, 'a',...
    params.plot_mx, params.plot_my, ...
    params.maux, get_kmax(params.space_order), components);
end

function num = get_aux_num(name)
    num = sscanf(name,'aux%d');
    if(~isnumeric(num)) num=-1; end;
end

% obtain a list of names that need to be defined
% so that the eval_string can be evaluated
%
% eval_string might be something like '-Pi13x-Pi23y' '-Pi23y-Pi13x'
function [var_names, first_name] = extract_var_names(c, eval_string)
  all_names = extract_names(c, eval_string);
  %
  %aux_names = {};
  %for i=1:numel(all_names)
  %  if(strncmp(all_names{i},'aux',3));
  %    num = get_aux_num(all_names{i});
  %    if(num > 0) aux_names = [aux_names all_names{i}]; end;
  %  end
  %end
  %
  % remove names of functions
  %special_names = {'log','exp','tanh','cosh','sinh'};
  %names=setdiff(names,special_names);
  %
  % should handle generic names like q1 or a1 specially
  %
  % restrict to valid names
  var_names=intersect(all_names,c.s.varMap.name_arr);
  if(numel(var_names)==0)
    error(['no valid name found among: ' all_names{:}]);
  end
  first_name = var_names{1};
end

function display_string = replace_var_names(c,var_names,eval_string)
  % construct array of replacement strings for var_names
  %
  rep_strings = cell(size(var_names));
  display_names = c.s.plot_params.display_names;
  varIdx_s = c.s.varMap.varIdx_s;
  for i=1:numel(var_names)
    var_name = var_names{i};
    rep_string = display_names{varIdx_s.(var_name)};
    rep_strings{i} = ['(' rep_string ')'];
  end

  % replace each string in var_names with its equivalent in rep_strings
  display_string = regexprep(eval_string,var_names,rep_strings);
end
