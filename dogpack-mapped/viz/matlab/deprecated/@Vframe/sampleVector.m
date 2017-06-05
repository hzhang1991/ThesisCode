
% components: an index, array of indices, component name, or cell array of names
% point_type = 1: gauss quadrature points (implies cell_pt_ratio = space_order)
% point_type = 2: cell_pt_ratio = pts_per_cell (default)
% point_type = 3: cell_pt_ratio = cells_per_pt (default:  space_order)
function values = sampleVector(frame, components, point_type, cell_pt_ratio)
  % if components is a cell array convert it to an array of indices

  components=get_componentIndices(frame,components);
  componentNames = get_componentNames(frame);
  if(~exist('point_type')) point_type = 2; end
  simulation = get_simulation(frame);
  space_order = get_space_order(simulation);
  if(~exist('cell_pt_ratio'))
    cell_pt_ratio = space_order;
    cell_pt_ratio = space_order;
  end
  % 3 is not yet supported
  assert(any(point_type==1:2));

  vector_components = frame.s.vector(:,:,:,components);
  values = sample_state2(vector_components, space_order, point_type, cell_pt_ratio);
end
