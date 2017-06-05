% constructor for Vframe ("vector frame") class
%
% num_frame: vframe number
%
% simulation: either the name of an output directory
%   or an instance of Parameters
%
function vframe=Vframe(num_frame, simulation, vectorName, components)

  % make sure that simulation has the data we need
  %
  if(nargin<4)
    vectorName='q';
  end
  assert(vectorName=='q'||vectorName=='a');
  %
  if(~exist('simulation'))
    global simulation;
    if(isempty(simulation))
      setoutputdir('output');
    end
  else
    if(ischar(simulation))
      simulation=Simulation(simulation);
    end
  end
  %
  assert(isa(simulation,'Parameters'));

  parameters = get_parameters(simulation);
  outputdir    = get_outputdir   (parameters);
  datafmt      = get_datafmt     (parameters);
  mx           = get_mx          (parameters);
  my           = get_my          (parameters);
  space_order  = get_space_order (parameters);
  if(nargin<3)
    if(vectorName=='a')
      meqn = get_maux(parameters);
      components=1:maux;
    else
      meqn = get_meqn(parameters);
      components=1:meqn;
    end
  end
  if(vectorName=='a')
    all_componentNames = get_a_componentNames(parameters);
  else
    all_componentNames = get_q_componentNames(parameters);
  end

  vector=[];
  if(numel(components)>0)
    %basefilename = [outputdir '/' vectorName sprintf('%04d',num_frame)];
    %[vector,time]=read_state2_cart(datafmt, basefilename, vectorName, ...
    %  mx, my, meqn, get_kmax(space_order), components);
    [vector,time]=read_state2_cart(datafmt, outputdir, num_frame, vectorName, ...
      mx, my, meqn, get_kmax(space_order), components);
  end
  %
  s = struct;
  s.componentNames = SubComponentNames(all_componentNames, components);
  s.vector=vector;
  s.time=time;
  s.parameters=parameters;

  vframe = struct;
  vframe.s = s;
  %
  vframe=class(vframe,'Vframe');
end
