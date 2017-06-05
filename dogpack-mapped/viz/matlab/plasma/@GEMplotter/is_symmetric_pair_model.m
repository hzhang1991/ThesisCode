
function bool = is_symmetric_pair_model(c)
  model_name = c.s.params.model_name;
  bool = strcmp(model_name, 'p10') || strcmp(model_name, 'p05');
end
