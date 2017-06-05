% The user can call this function while in a dogpack directory
% to ensure that all 
%
% should rename this something like setdogpath.m
% so it won't conflict with other matlab environments.
%
% could also create a top-level 'setdog' script which
% initializes the dogpack environment
% (e.g. by calling setdogpath and setoutputdir).
%
function setpath(varargin)
  % set dogpack matlab directories
  %[status, DOGPACK] = system('echo -n $DOGPACK');
  setdogpath;
  setdogdimpath(varargin{:});
end
