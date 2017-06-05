% I should change this to look at the $DOGPACK/config/*_defaults.ini
% files to get the default values
%
function params = get_GEM_params(outputdir)
  params = get_dogpack_params(outputdir); % method
  % plasma model parameters are GEM parameters that
  % make sense independent of the GEM problem
  params = add_plasma_model_params(outputdir, params); % model
  params = add_GEM_params(outputdir, params); % problem
end

function params = add_plasma_model_params(outputdir, params)
  parameters_ini = [outputdir '/parameters.ini'];
  param_names_and_defaults = { ...
    ... % the model name may also be a GEM method parameter
    'model_name', 'none', ...
    ... % the light speed and masses are both model and problem parameters
    'cs_light', -1, ...
    'gamma', 5./3., ...
    ... % masses
    'mass_ratio', -1, ...
    'spc_mass_mode', 'total_mass', ...
    'ion_mass', -1, ...
    'total_mass',  1, ...
    ... % isotropization
    'iso_period_type', 'default', ...
    'ion_iso_period', -1, ...
    'elc_iso_period', -1, ...
    'base_iso_period', -1, ...
    ... % resistivity
    'slowing_period', -1, ...
    'trigger_mui3', -1, ...
    'slowing_rate_slope', -1};
  param_names=param_names_and_defaults(1:2:end);
  param_defaults=param_names_and_defaults(2:2:end);
  for(i=1:length(param_names))
    the_name = param_names{i};
    the_default = param_defaults{i};
    params.(the_name) = param_defaults{i};
  end
  params = read_more_params_from_ini(parameters_ini, params, param_names{:});
%    ... % the model name may also be a GEM method parameter
%    'model_name', ...
%    ... % the light speed and masses are both model and problem parameters
%    'cs_light', ...
%    'gamma' , ...
%    'mass_ratio','spc_mass_mode','ion_mass','total_mass',   ... % masses
%    'ion_iso_period', 'elc_iso_period',                     ... % isotropization
%    'slowing_period', 'trigger_mui3', 'slowing_rate_slope'  ... % resistivity
%    );
    %
  params.outputdir = outputdir;

  % set derived parameters
  %
  params.one_over_epsilon = params.cs_light*params.cs_light;
  if(strcmp(params.spc_mass_mode,'total_mass'))
    if(params.total_mass == -1) params.total_mass=1.; end
    params.elc_mass = params.total_mass/(1.+params.mass_ratio);
    params.ion_mass = params.mass_ratio*params.elc_mass;
  elseif(strcmp(params.spc_mass_mode,'ion_mass'))
    if(params.ion_mass == -1) params.ion_mass=1.; end
    params.elc_mass = params.ion_mass/params.mass_ratio;
    params.total_mass = params.ion_mass+params.elc_mass;
  else
    error(['invalid value of params.spc_mass_mode:' params.spc_mass_mode]);
  end
  %
  % set spc_base_iso_period
  if(strcmp(params.iso_period_type,'default'))
    if(params.base_iso_period > 0)
      params.iso_period_type = 'det';
    else
      params.iso_period_type = 'constant';
    end
  end
  switch params.iso_period_type
  case {'det','trace'}
    params.ion_base_iso_period = params.base_iso_period;
    params.elc_base_iso_period = params.base_iso_period;
  case 'constant'
  otherwise
    error(['unsupported iso_period_type: ' species]);
  end
  %
  % set model name
  %
  if(strcmp(params.model_name,'none')) % do nothing
    params.model_string = 'plasma model';
  elseif(strcmp(params.model_name,'g10'))
    params.model_string = '10-moment two-fluid plasma';
  elseif(strcmp(params.model_name,'g05'))
    params.model_string = '5-moment two-fluid plasma';
  elseif(strcmp(params.model_name,'p10'))
    params.model_string = '10-moment symmetric pair plasma';
  elseif(strcmp(params.model_name,'p05'))
    params.model_string = '5-moment symmetric pair plasma';
  elseif(strcmp(params.model_name,'i10e5'))
    params.model_string = '10-moment ions, 5-moment electrons';
  elseif(strcmp(params.model_name,'mhd'))
    params.model_string = 'MHD';
  else
    error(['unsupported model_name:' model_name]);
  end
  
end

function params = add_GEM_params(outputdir, params);
  parameters_ini = [outputdir '/parameters.ini'];
  params = read_more_params_from_ini(parameters_ini, params, ...
    ... % GEM method parameters
    'enforce_symmetry', 'mesh_is_virtual', ...
    'enforced_symmetry', ... % deprecated
    ... % GEM problem parameters
    'BCs', ...
    'B_0', 'n_0', 'domain_scaling', 'sheet_thickness',  ... % GEM scaling
    'B_guide', ...
    'temp_ratio'                                        ... % temp ratio
    );
  out_parameters_ini = [outputdir '/out_parameters.ini'];
  params = read_more_params_from_ini(out_parameters_ini, params, ...
    'neg_x_is_virtual', ...
    'neg_y_is_virtual', ...
    'enforcing_rotational_symmetry');

  if(~isfield(params,'B_guide'))
    params.B_guide = 0.;
  end
  if(~isfield(params,'n_0'))
    params.n_0 = 1.;
  end
  if(~isfield(params,'sheet_thickness'))
    params.sheet_thickness = 1.;
  end
  if(~isfield(params,'enforced_symmetry'))
    params.enforced_symmetry = params.neg_x_is_virtual ...
                           + 2*params.neg_y_is_virtual;
  else
    warning('using deprecated parameter enforced_symmetry');
    params.neg_x_is_virtual = bitand(params.enforced_symmetry,1);
    params.neg_y_is_virtual = ~~bitand(params.enforced_symmetry,2);
    params.enforcing_rotational_symmetry = 0;
  end

  % for backwards compatibility with data
  if(params.out_ini==0)
    % override dogpack parameters based on GEM parameters
    %
    params.xhigh=4*pi*params.domain_scaling;
    params.xlow=-params.xhigh;
    params.yhigh=2*pi*params.domain_scaling;
    params.ylow=-params.yhigh;
    if(bitand(params.enforced_symmetry,1)==1)
      params.xlow=0;
    end
    if(bitand(params.enforced_symmetry,2)==2)
      params.ylow=0;
    end
    % bit 4 means mx and my describe actual computational mesh
    if(bitand(params.enforced_symmetry,5)==1)
      assert(bitand(params.mx,1)==0); % mx is even
      params.mx = params.mx/2;
    end
    if(bitand(params.enforced_symmetry,6)==2)
      assert(bitand(params.my,1)==0); % my is even
      params.my = params.my/2;
    end
    % params = set_derived_dogpack_params(params);
    params.dx = (params.xhigh - params.xlow)/params.mx;
    params.dy = (params.yhigh - params.ylow)/params.my;
    params.plot_dx = params.dx;
    params.plot_dy = params.dy;
  end

  % set nondimensionalizing quantities
  %
  params.alfven_speed = params.B_0/sqrt(params.n_0*params.total_mass);
  % params.ion_alfven_speed = params.B_0/sqrt(params.n_0*params.ion_mass);
  params.ion_gyrofreq = params.B_0/params.ion_mass;
  % gyroradius is thermal velocity divided by angular gyrofreq
  % thermal velocity is 
  % for the GEM problem typical pressure is (half) B_0^2,
  % so a typical gyroradius is (half) the ion skin dept
  if(strcmp(params.model_name, 'mhd'))
    params.ion_skindepth = 1;
  else
    params.ion_skindepth = sqrt(params.ion_mass/params.n_0);
  end
  % this is the same as params.scale.E_0 only if total_mass is chosen below.
  params.E_0 = params.B_0*params.alfven_speed;

  % modify commenting to choose Alfven speed or ion Alfven speed
  % for nondimensionalization of plotted output.
  % (note that the scale of the domain is determined by
  % the solver and is not modified when displaying data.)
  %
  scale = make_scale(params,params.total_mass);
  %scale = make_scale(params,params.ion_mass);

  params.scale=scale;

  % check parameters
  %
  if(~strcmp(params.model_name, 'mhd'))
    assert(test_equal(params.total_mass, params.elc_mass+params.ion_mass));
    assert(test_equal(params.ion_mass/params.elc_mass,params.mass_ratio));
    assert(test_equal(params.mass_ratio,params.temp_ratio^2));
    % if(params.ion_iso_period<0) assert(params.elc_iso_period<0);
    % else
    %   assert(test_equal(params.ion_iso_period, ...
    %                     params.elc_iso_period*params.temp_ratio));
    % end
  end
end

% The original GEM paper [Birn01] says that the normalization
% of the space scales and timescales of the system is chosen
% to be the ion inertial length (a.k.a. skin depth) and the
% ion cyclotron frequency (a.k.a. gyrofrequency), where the
% inertial length is evaluated with the density n_0 and
% the ion gyrofrequency is evaluated at the peak magnetic
% field. They then say that in these units the velocities
% are normalized to the Alfven speed, which actually is
% not quite exactly true because the electrons have some
% mass. Specifically, the formula for the Alfven speed is
% proportional to B_0/sqrt(rho), and rho = rho_i + rho_e.
% (I say "proportional to" to evade the factor of 1=sqrt(\mu_0)
% or 1/sqrt(4*\pi) depending on SI or gaussian units.)
% For the case of pair plasma this difference becomes more
% pronounced.
% 
% Following the precedent of Bessho and Bhattacharjee, whose
% results we seek to replicate, we use the ion mass (and hence
% the "ion Alfven speed") to nondimensionalize space and time,
% but when reporting output values we instead
% nondimensionalize by total mass.
%
function scale = make_scale(params,m_0)
  % The consistent scales that would be used to rescale the input:
  %
  scale.m_0 = m_0;
  scale.n_0 = params.n_0;
  %scale.x_0 = params.ion_skindepth;
  scale.B_0 = params.B_0;
  %
  % derived scales
  %
  scale.u_0 = scale.B_0/sqrt(scale.n_0*scale.m_0);
  scale.t_0 = scale.m_0/scale.B_0; % "gyroperiod"
  scale.x_0 = scale.u_0*scale.t_0;
  %scale.u_0 = scale.x_0/scale.t_0; % "Alfven speed"
  scale.rho_0 = scale.n_0*scale.m_0;
  scale.M_0 = scale.rho_0*scale.u_0;
  scale.P_0 = scale.M_0*scale.u_0;
  scale.E_0 = scale.B_0*scale.u_0;
  scale.J_0 = scale.n_0*scale.u_0;
  scale.T_0 = scale.P_0/scale.n_0;
end

