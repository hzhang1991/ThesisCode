
function plotScalar(frame, component)
  componentNames = get_componentNames(frame);
  if(ischar(component))
    component=get_index(componentNames, component);
  end
  assert(numel(component)==1);
  name = get_name(componentNames, component);
  parameters = get_parameters(frame);
  assert(get_ndims(parameters)==2);

  % access component
  vector = get_vector(frame);
  vector_component = vector(:,:,:,component);

  % sample values
  space_order = get_space_order(parameters);
  pts_per_cell = space_order;
  values = sample_state2(vector_component, space_order, 1, pts_per_cell);

  % plot values
  grid = get_grid(parameters);
  axis_array = [...
    grid.xlow ,...
    grid.xhigh,...
    grid.ylow ,...
    grid.yhigh];
  mxpts = grid.mx*pts_per_cell;
  mypts = grid.my*pts_per_cell;
  xl = grid.xlow+(0:mxpts)*(grid.dx/pts_per_cell);
  yl = grid.ylow+(0:mypts)*(grid.dy/pts_per_cell);
  fig = component;
  plot_scalar2(fig, values, xl, yl, get_time(frame), axis_array, name);
end
