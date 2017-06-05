function parameters = read_params_from_ini(ini_fileName, varargin)
%=========================================================================
% Define variables from .ini file ignoring section labels
% ini_fileName: ini file (including .ini extension)
%==========================================================================
  parameters = struct;
  parameters = read_more_params_from_ini(ini_fileName, parameters, varargin{:});
end
