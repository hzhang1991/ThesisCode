% constructor for Simulation class
% This class is a broker to access simulation output.
%
% for now it just inherits from and thus is identical to
% the Parameters class.
%
function simulation=Simulation(varargin)
  %
  simulation=struct;
  % pass all arguments through to the parent constructor
  parameters = Parameters(varargin{:});
  %simulation.grid = Grid(get_cartesian(parameters));
  % inherit from parameters
  simulation=class(simulation,'Simulation',parameters);
end
