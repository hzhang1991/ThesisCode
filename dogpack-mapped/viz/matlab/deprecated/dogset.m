
function dogset(varargin)
  %profile on;
  % set path components that don't depend on parameters.ini
  dogsetpath
  dogsetoutputdir(varargin{:});
  addpath('.', './matlab')
  %profile off;
  %profile viewer;
end

% set the default output directory
% and associated information
%
function dogsetoutputdir(outputdir_in)
  setoutputdir(outputdir_in);
  if(nargin<1)
    outputdir_in='output'
  end
  global simulation;
  global outputdir;
  outputdir=outputdir_in;
  simulation = Simulation(outputdir);
  
  % reset components of MATLABPATH that could depend on 
  % the simulation parameters
  
  matlab_addpath = get_matlab_addpath(simulation);
  if(~isempty(matlab_addpath))
    eval(matlab_addpath);
  end
  dogset_dim_path(get_ndims(simulation));
end

function dogset_dim_path(dim)
  % if dim not provided, try to infer from pwd
  % (or could obtain from parameters.ini).
  % ...
  if(nargin<1)
    dim=get_ndims_space();
  end
  addpath([getenv('DOGPACK') '/matlab/' num2str(dim) 'd']);
end

function dim=get_ndims_space()
  global simulation;
  if(isa(simulation,'Simulation'))
    dim=get_ndims(simulation);
  else
    dim=1; % the default
    parsed_path = textscan(pwd, '%s', 'delimiter', '/');
    for i=1:numel(parsed_path{1})
      if(strcmp(parsed_path{1}(i),'2d'))
        dim=2;
      elseif(strcmp(parsed_path{1}(i),'1d'))
        dim=1;
      elseif(strcmp(parsed_path{1}(i),'3d'))
        dim=3;
      end
    end
  end
end
