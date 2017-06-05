function plotall(outputdir_in)
  setarg_outputdir;
  plot_recon(outputdir);
  model_name = read_param_from_ini([outputdir_in '/parameters.ini'],...
      'model_name');
  plot_xpoint(outputdir,40,1,'i');
  command=['plotout' model_name '(outputdir);'];
  disp(['evaluating command: ' command]);
  eval(command);
end

% This function moved from viz/matlab/read_param_from_ini because plotall.m
% was the only item calling this. Oct. 3, 2013 (-DS).
%
function optionVal = read_param_from_ini(ini_fileName, optionName)
% could accelerate this for a scalar parameter by making it return
% immediately upon finding match

  parameters = read_params_from_ini(ini_fileName, optionName);
  if(~isfield(parameters,optionName))
    error(['no option named ' optionName ' in file ' ini_fileName]);
  end
  optionVal = parameters.(optionName);
end
