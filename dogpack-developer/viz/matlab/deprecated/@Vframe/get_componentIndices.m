
function components = get_componentIndices(frame,components)
  if(ischar(components))
    components=get_index(frame.componentNames, components);
  elseif(iscell(components))
    nameArray=components;
    components=zeros(1,numel(nameArray));
    for i=1:numel(nameArray)
      % convert name to index
      components(i)=get_index(frame.componentNames,nameArray{i});
    end
  end
end

