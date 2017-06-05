% constructor for Parameters class
%
function v=Parameters(outputdir)
  %
  assert(isdir(outputdir),['no such directory: ' outputdir]);
  %
  % fetch the parameters
  %
  % should extract only the parameters we actually
  % want to use for plotting?
  parameter_file = [outputdir '/parameters.ini'];
  parameters = read_ini(parameter_file);
  %
  % translate old names to new names
  %
  assert(isfield(parameters, 'dogParams'),...
    ['missing required section [dogParams] in parameter file ' parameter_file]);
  % for backwards compatibility
  if(isfield(parameters.dogParams,'method'))
    method = parameters.dogParams.method;
    parameters.dogParams.space_order = method(1);
    parameters.dogParams.maux = method(6);
    rmfield(parameters.dogParams, 'method');
  end
  %
  % add parameter information
  %
  parameters.outputdir=outputdir;
  %
  if(strcmp(parameters.dogParams.mesh_type,'Cartesian'))
    assert(isfield(parameters,'grid'));
    grid = parameters.grid;
    space_order = parameters.dogParams.space_order;
    ndims = parameters.dogParams.ndims;
    grid.kmax = get_kmax(space_order,ndims);
    grid.dx = (grid.xhigh-grid.xlow)/grid.mx;
    grid.dy = (grid.yhigh-grid.ylow)/grid.my;
    parameters.grid = grid;
  end
  %
  % create q_names
  %
  meqn = parameters.dogParams.meqn;
  assert(isfield(parameters, 'componentParams'));
  assert(isfield(parameters.componentParams, 'q_names'));
  if(isfield(parameters, 'componentParams')...
  && isfield(parameters.componentParams, 'q_names'))
    assert(numel(parameters.componentParams.q_names) == meqn);
    q_names = parameters.componentParams.q_names;
  else
    % use default q_names: q1, q2, etc.
    %parameters.componentParams=struct;
    q_names=cell(1,meqn);
    for i=1:meqn
      q_names{i}=['q' num2str(i)];
    end
  end
  %
  % create a_names
  %
  maux = parameters.dogParams.maux;
  if(isfield(parameters, 'componentParams')...
  && isfield(parameters.componentParams, 'a_names'))
    assert(numel(parameters.componentParams.a_names) == maux);
    a_names = parameters.componentParams.a_names;
  else
    % use default a_names: a1, a2, etc.
    %parameters.componentParams=struct;
    a_names=cell(1,maux);
    for i=1:maux
      a_names{i}=['a' num2str(i)];
    end
  end

  % create mapping from names to components
  parameters.q_componentNames = ComponentNames(q_names);
  parameters.a_componentNames = ComponentNames(a_names);
  v=struct;
  v.parameters=parameters;
  v=class(v,'Parameters');
end

function name2idx = get_name2idx(names);
  meqn = numel(names);
  name2idx = struct;
  for i=1:meqn
    name2idx.(names{i})=i;
  end
end
