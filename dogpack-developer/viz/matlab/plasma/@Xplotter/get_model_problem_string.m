
function model_problem_string = get_model_problem_string(xplotter, params,do_rescale)
    model_problem_string = [ ...
      get_model_string(xplotter, params)  ...
      ', ' get_full_problem_string(xplotter, params,do_rescale)];
end

