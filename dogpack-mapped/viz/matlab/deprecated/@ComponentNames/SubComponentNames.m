function c_out = SubComponentNames(c_in, components)
  % if components is cell array convert it to indices
  components=get_indices(c_in,components);
  nameArray = get_nameArray(c_in);
  c_out = ComponentNames(nameArray(components));
end
