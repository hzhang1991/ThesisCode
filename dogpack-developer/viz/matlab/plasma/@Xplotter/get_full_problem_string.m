
function full_problem_string = get_full_problem_string(xplotter, params,do_rescale)
  if(do_rescale)
    full_problem_string = [ ...
      'm_i/m_e = ' num2str(params.mass_ratio) ...
      ', c/v_A = ', sprintf('%2.1f',params.cs_light/params.alfven_speed), ...
      ', ' get_rescale_problem_string(params) ...
    ' ' ];
  else
    full_problem_string = [ ...
      get_mass_string(params)  ... % comment out?
      ', ' get_problem_string(params) ... % comment out?
    ' ' ];
  end
end

function mass_str = get_mass_string(params)
  mass_ratio_str = ['m_i/m_e = ' num2str(params.mass_ratio)];
  % get mass_scale_str
  if(isfield(params,'spc_mass_mode') & ...
     strcmp(params.spc_mass_mode,'ion_mass'))
    ion_mass = 1.;
    if(isfield(params,'ion_mass')) ion_mass = params.ion_mass; end;
    mass_scale_str = ['m_i=' num2str(ion_mass)];
  else
    total_mass = 1.;
    if(isfield(params,'total_mass')) total_mass = params.total_mass; end;
    mass_scale_str = ['m_i+m_e=' num2str(total_mass)];
  end
  mass_str = [mass_ratio_str ', ' mass_scale_str];
end

function problem_str = get_problem_string(params)
  problem_str = [ ...
    'T_i/T_e = ' num2str(params.temp_ratio) ...
    'domain scaling = ' num2str(params.domain_scaling) ...
    ];
end

function problem_str = get_rescale_problem_string(params)
  nondim_Ly = params.domain_scaling/params.ion_skindepth;
  nondim_Lx = 2*nondim_Ly;
  nondim_sheet_thickness = params.sheet_thickness/params.ion_skindepth;
  problem_str = [ ...
    ... %'T_i/T_e = ' num2str(params.temp_ratio) ...
    ', 2\lambda = ' sprintf('%.2f',nondim_sheet_thickness) '\delta_i' ...
    ', L_x/(8\pi) = ' sprintf('%.2f',nondim_Ly) '\delta_i' ...
    ];
end

