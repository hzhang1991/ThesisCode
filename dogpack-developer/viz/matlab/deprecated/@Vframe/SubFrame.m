% select a subset of the existing components
function vframe = SubFrame(vframe, components)
  % if components is a cell array convert it to an array of indices
  if(iscell(components))
    components=get_indices(vframe.s.componentNames,components);
  end
  vframe.s.vector = vframe.s.vector(:,:,:,components);
  vframe.s.componentNames = SubComponentNames(vframe.s.componentNames, components);
end
