
function out = get_matlab_addpath(v)
  if(isfield(v.parameters.dogParams, 'matlab_addpath'))
    out = v.parameters.dogParams.matlab_addpath;
  else
    out = '';
  end
end

