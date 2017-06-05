
function out = is_symmetric_pair_model(c)
  model_name = c.s.params.model_name;
  out = strcmp(model_name, 'p20') || strcmp(model_name, 'p10') || strcmp(model_name, 'p05');
end
