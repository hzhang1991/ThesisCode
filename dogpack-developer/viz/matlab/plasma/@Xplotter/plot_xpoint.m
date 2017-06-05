
% species: 'i' or 'e' (ion or electron)
%function plot_xpoint(xplotter,outputdir_in, plots, final_time, species, ...
%    do_rescale, show_islands, figure_shift)
function plot_xpoint(xplotter,outputdir_in, varargin)

  setarg_outputdir;

  % set default values of keyword arguments
  %
  do_print=0;
  plots=1;
  species='i';
  figure_shift=100;
  do_rescale=3;
  show_islands=3;
  ylims=[-20 80];
  %
  % parse keyword arguments
  %
  idx=1;
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      %case 'dir'; idx=idx+1;
      %  outputdir_in=varargin{idx}; idx=idx+1;
      case {'print','do_print'}; idx=idx+1;
        do_print=1;
      case 'plots'; idx=idx+1;
        plots=varargin{idx};
        idx=idx+1;
      case 'species'; idx=idx+1;
        species=varargin{idx};
        idx=idx+1;
      case 'figure_shift'; idx=idx+1;
        figure_shift=varargin{idx}; idx=idx+1;
      case 'do_rescale'; idx=idx+1;
        do_rescale=varargin{idx}; idx=idx+1;
      case 'show_islands'; idx=idx+1;
        show_islands=varargin{idx}; idx=idx+1;
      case 'final_time'; idx=idx+1;
        final_time=varargin{idx}; idx=idx+1;
      case 'ylims'; idx=idx+1;
        ylims=varargin{idx}; idx=idx+1;
      otherwise
        error(['invalid keyword: ' varargin{idx}]);
    end
  end
  %
  % set derived values from bits of arguments
  %
  plot_struct.plot_accum = bitand(plots,1);
  plot_struct.plot_ohm = bitand(plots,2);
  plot_struct.plot_resistivity = bitand(plots,4);
  % first bit means rescale space-time.
  plot_struct.do_rescale_times = bitand(do_rescale,1);
  plot_struct.do_rescale_lengths = bitand(do_rescale,1);
  plot_struct.do_rescale_fluxes = bitand(do_rescale,2);
  plot_struct.use_gauss = bitand(do_rescale,4);
  plot_struct.show_percent_flux = 0;
  plot_struct.show_islands = 0;
  plot_struct.ylims = ylims;
  if(show_islands)
    fluxes = get_fluxes(Xplotter(),outputdir);
  else
    fluxes.show_central_island_flux = 0;
    fluxes.show_other_island_flux = 0;
  end

  % check arguments
  assert(strcmp(species,'i')|strcmp(species,'e')|strcmp(species,'t'));

  %params = get_GEM_params(outputdir);
  [outputdir_shouldname, params] = make_outputdirname(outputdir);

  [t, rho_i, M3i,  P13i_x, P23i_y, E3, ...
   mu3i, inertial_term, divPs_term, Ps13_x_term, Ps23_y_term, ...
      resistive_term, residual] = ...
    get_ohms_law_terms(Xplotter(), species, params, outputdir);

  %%%%%%%%%%%%%%%%%%%%
  %%% plot data
  %%%%%%%%%%%%%%%%%%%%

  %fonts.ticlabelfontsize = 30;
  %fonts.annotationfontsize = 34;
  %fonts.xlabelfontsize = 32;
  %fonts.ylabelfontsize = 32;
  %fonts.titlefontsize = 36;
  %fonts.legendfontsize = 40;

  fonts.ticlabelfontsize = 14;
  fonts.annotationfontsize = 17;
  fonts.xlabelfontsize = 18;
  fonts.ylabelfontsize = 16;
  fonts.titlefontsize = 18;
  fonts.legendfontsize = 16;

  if(~exist('final_time','var'))
    final_time = t(end);
  end

  plot_struct.outputdir_shouldname = outputdir_shouldname;
  plot_struct.fonts = fonts;
  plot_struct.figure_shift = figure_shift;
  plot_struct.final_time = final_time;

  if(plot_struct.do_rescale_times)
    plot_struct.final_time = plot_struct.final_time*params.ion_gyrofreq;
  end
  using_resistive_term = get_using_resistive_term(params);

  if(plot_struct.plot_accum)
    plot_accum(plot_struct, t, E3, divPs_term, Ps13_x_term, Ps23_y_term, ...
      inertial_term, mu3i, ...
      using_resistive_term, resistive_term, residual, fluxes, params, species);
    if(do_print); print_plot(figure_shift, params,1); end;
  end
  if(plot_struct.plot_ohm)
    plot_ohm(plot_struct, t, E3, divPs_term, Ps13_x_term, Ps23_y_term, ...
      inertial_term, ...
      using_resistive_term, resistive_term, residual, params, species);
    if(do_print); print_plot(figure_shift, params,2); end;
  end
  if(plot_struct.plot_resistivity)
    residual_slowing_rate = residual./mu3i;
    %effective_resistivity = resistive_term./u3i;
    plot_resistivity(plot_struct, t, residual_slowing_rate, ...
      mu3i, M3i, params, species);
  end
  
end
  
function [t, rho_s, M3s,  P13s_x, P23s_y, E3, ...
    mu3s, inertial_term, divPs_term, Ps13_x_term, Ps23_y_term, ...
    resistive_term, residual] = ...
get_ohms_law_terms(Xplotter, species, params, outputdir)

  % set species masses
  %[ion_mass, elc_mass] = get_species_masses(params);
  total = 0;
  switch species
    case 'i'
      spc_mass = params.ion_mass;
    case 'e'
      spc_mass = params.elc_mass;
    case 't' % total
      spc_mass = params.ion_mass + params.elc_mass;
      total = 1;
    otherwise
      error(['invalid species: ' species]);
  end

  if(total==0)
    % read xpoint data
    [t, rho_s, M3s, P13s_x, P23s_y, E3] = read_xpoint_data(Xplotter(),...
      outputdir, species);
    % calculate the ohm's law terms
    [mu3s, inertial_term, divPs_term, Ps13_x_term, Ps23_y_term, ...
        resistive_term, residual] = ...
      calculate_ohms_law_terms(Xplotter(), t, rho_s, M3s,  P13s_x, P23s_y, E3, ...
        species, spc_mass, params);
  else
    % for combined ohm's law just average the two "ohm's laws"
    %
    % read xpoint data
    [ti, rho_i, M3i, P13i_x, P23i_y, E3i] = read_xpoint_data(Xplotter(),...
      outputdir, 'i');
    [t, rho_e, M3e, P13e_x, P23e_y, E3] = read_xpoint_data(Xplotter(),...
      outputdir, 'e');
    assert(all(ti==t));
    assert(all(E3i==E3));
    % calculate the ohm's law terms
    [mu3i, inertial_term_i, divPi_term, Pi13_x_term, Pi23_y_term, ...
        resistive_term_i, residual_i] = ...
      calculate_ohms_law_terms(Xplotter(), t, rho_i, M3i,  P13i_x, P23i_y, E3, ...
        'i', params.ion_mass, params);
    [mu3e, inertial_term_e, divPe_term, Pe13_x_term, Pe23_y_term, ...
        resistive_term_e, residual_e] = ...
      calculate_ohms_law_terms(Xplotter(), t, rho_e, M3e,  P13e_x, P23e_y, E3, ...
        'e', params.elc_mass, params);

    % average the output quantities from calculate_ohms_law_terms
    mu3s             = 0.5.*(mu3i             + mu3e);
    inertial_term    = 0.5.*(inertial_term_i  + inertial_term_e);
    divPs_term       = 0.5.*(divPi_term       + divPe_term);
    Ps13_x_term      = 0.5.*(Pi13_x_term      + Pe13_x_term);
    Ps23_y_term      = 0.5.*(Pi23_y_term      + Pe23_y_term);
    resistive_term   = 0.5.*(resistive_term_i + resistive_term_e);
    residual         = 0.5.*(residual_i       + residual_e);

    % sum other quantities
    rho_s = rho_i + rho_e;
    M3s = M3i + M3e;
    P13s_x = P13i_x + P13e_x;
    P23s_y = P23i_y + P23e_y;

    % % combine the corresponding quantities
    % % to calculate the full ohm's law terms
    % %
    % rho_s = rho_i+rho_e;
    % M3s = M3i + M3e;
    % P13s_x = P13i_x + P13e_x;
    % P23s_y = P23i_y + P23e_y;
    % Ps13_x_term = (P13i_x./rho_i).*params.ion_mass -  ...
    %               (P13e_x./rho_e).*params.elc_mass;
    % Ps23_y_term = (P23i_y./rho_i).*params.ion_mass -  ...
    %               (P23e_y./rho_e).*params.elc_mass;
    % divPs_term = (Ps13_x_term + Ps23_y_term)./2.;
    % ni = rho_i./params.ion_mass;
    % ne = rho_e./params.elc_mass;
    % mu3s = (M3i./ni - M3e./ne)./2.;
  end

end

function print_plot(figure_shift, params, fig)
  figure(fig+figure_shift);

  if(fig==1)
    fig_str = 'accum';
  elseif(fig==2)
    fig_str = 'ohm';
  else
    error(['invalid value of figure: ' num2str(figure)]);
  end

  % make iso_period_str
  %
  iso_period_str='';
  if(strcmp(params.model_name,'p10') || strcmp(params.model_name,'p20'))
    iso_period = num2str(params.ion_iso_period);
    if(params.ion_iso_period < 0)
      iso_period = 'infty';
    end
    % replace periods with underscores so latex doesn't get confused.
    %iso_period = regexprep(iso_period, '\.', '_');
    iso_period_str = ['_iso_period=' iso_period];
  end
  cfl_str = ['.cfl=' num2str(params.cflv(2))];
  
    basename = [getenv('HOME') '/figures'];
    filename = [basename '/' fig_str '_' ...
      params.model_name '_' num2str(params.mx) 'x' num2str(params.my) ...
      cfl_str  iso_period_str];
    filename = regexprep(filename, '\.','_');
    disp(['print to file: ' filename]);
    print('-depsc2',filename);
end

function qminus=get_qminus(species)
  qminus='';
  if(species=='i')
    qminus='-';
  end
end
  
function out=integrate(t,f)
  fh = (f(2:end) + f(1:end-1))*0.5;
  dt = t(2:end)-t(1:end-1);
  out=zeros(size(f));
  out(2:end) = cumsum(fh.*dt);
end

function using_resistive_term = get_using_resistive_term(params)
  using_resistive_term = 0;
  if(params.slowing_period > 0)
    using_resistive_term = 1;
  elseif(params.slowing_rate_slope > 0)
    using_resistive_term = 1;
  end
end

% need to fix this to handle rescaling
function slowing_string = get_slowing_string(params, ...
  do_rescale_times,do_rescale_fluxes)

  assert(params.slowing_period ~= 0.);
  slowing_string = '';
  if(params.slowing_period > 0)
    %if(do_rescale_times)
    %  slowing_period = params.slowing_period*params.ion_gyrofreq;
    %else
      slowing_period = params.slowing_period;
    %end
    slowing_string = ['slowing period: ' num2str(slowing_period)];
  elseif(params.slowing_rate_slope > 0)
    %if(do_rescale_fluxes)
    %  flux_scaling = (params.ion_gyrofreq / params.E_0);
    %  trigger_mui3 = params.trigger_mui3 * flux_scaling;
    %  slowing_rate_slope = params.slowing_rate_slope * flux_scaling;
    %else
      trigger_mui3 = params.trigger_mui3;
      slowing_rate_slope = params.slowing_rate_slope;
    %end
    slowing_string = ['trigger m_i u_3 = ' num2str(trigger_mui3) ...
      ', slowing rate slope = ' num2str(slowing_rate_slope)];
  end
end

function plot_ohm(plot_struct, t, E3, divPs_term, Ps13_x_term, Ps23_y_term, ...
  inertial_term, ...
  using_resistive_term, resistive_term, residual, params, species)

  fonts = plot_struct.fonts;

  % rescale fluxes if requested
  if(plot_struct.use_gauss)
    y_label = 'E_{gauss}';
  else
    y_label = 'E_{SI}';
  end
  if(plot_struct.do_rescale_fluxes)
    E_rescale_factor = 1./params.E_0;
    if(E_rescale_factor ~= 1.)
      E3 = E3*E_rescale_factor;
      divPs_term = divPs_term*E_rescale_factor;
      Ps13_x_term = E_rescale_factor*Ps13_x_term;
      Ps23_y_term = E_rescale_factor*Ps23_y_term;
      inertial_term = inertial_term*E_rescale_factor;
      resistive_term = resistive_term*E_rescale_factor;
      residual = residual*E_rescale_factor;
    end
    if(plot_struct.use_gauss)
      y_label = [ y_label 'c/(B_0 v_A)' ];
    else
      y_label = [ y_label '/(B_0 v_A)' ];
    end
  end
  % rescale time if requested
  if(plot_struct.do_rescale_times)
    t = t*params.ion_gyrofreq;
  end

  % plot Ohm's law terms at X-point
  %
  figure(2+plot_struct.figure_shift);
  clf;
  hold on;
  p1=plot(t,-E3             ,'b-'); set(p1,'linewidth',2);
  p1=plot(t,-divPs_term     ,'r-'); set(p1,'linewidth',3);
  p1=plot(t,-inertial_term  ,'k-'); set(p1,'linewidth',2);
  if(using_resistive_term)
  p1=plot(t,-resistive_term ,'g-'); set(p1,'linewidth',2);
  end
  p1=plot(t,-residual       ,'m-'); set(p1,'linewidth',1);
  p1=plot(t,-Ps13_x_term); set(p1,'linewidth',2,'color',[.6 .2 0],'linestyle',':');
  p1=plot(t,-Ps23_y_term); set(p1,'linewidth',2,'color',[.2 .6 0],'linestyle','-.');
  hold off;
  set(gca,'fontsize',fonts.ticlabelfontsize,'fontweight','bold');
  %t1=title([ '"Ohm''s law" terms at the X-point [DoGPack]']);
  t1=title([ 'Ohm''s law terms at the X-point [DoGPack]']);
  set(t1,'fontsize',fonts.titlefontsize,'fontweight','bold');
  xlim([t(1) plot_struct.final_time]);
  lo = -max([max(E3),max(divPs_term),max(inertial_term),max(resistive_term),max(residual)]);
  hi = -min([min(E3),min(divPs_term),min(inertial_term),min(resistive_term),min(residual)]);
  ylim([lo,hi]*1.1);
  %ylim([-.2 .2].*(params.domain_scaling*E_rescale_factor));
  axis on; box on; grid off;
  hold off
  %set(gca,'xtick',0:100:500);
  %set(gca,'ytick',-0.6:0.02:0.6);
  set(gca,'xtick',0:5:80);
  %ylim([-.16,.2]);
  grid on;
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  qminus=get_qminus(species);
  species_idx = ['_' species];
  if(species=='t')
    species_idx='';
  end
  labels = { ...
     [qminus 'E_z'], ...
     [qminus 'div(P_i-P_e)/(e(n_i+n_e))'], ...
     ... %[qminus 'div(P' species_idx ')_z/(en' species_idx ')'],...
     [qminus '(m_id_tu_i-m_ed_tu_e)/(2e)'],...
     ... %[qminus 'd_t(u' species_idx ')_z  m' species_idx '/e'],...
      };
  if(using_resistive_term)
    labels = [labels, 'resistive term', 'residual'];
  else
    labels = [labels, 'residual'];
  end
  labels = [labels, [qminus '(P_i-P_e)_x_z_,_x/(e(n_i+n_e))']];
  labels = [labels, [qminus '(P_i-P_e)_y_z_,_y/(e(n_i+n_e))']];
  %labels = [labels, ['(P' species_idx ')_x_z_,_x/(en' species_idx ')']];
  %labels = [labels, ['(P' species_idx ')_y_z_,_y/(en' species_idx ')']];
  fontsize = 13; %fonts.legendfontsize;
  l1 = legend(labels,...
    'location','northwest', ...
    'fontsize',fontsize);

  xlabel(get_time_label(Xplotter(),plot_struct.do_rescale_times),...
     'fontsize',fonts.xlabelfontsize,'fontweight','bold');
  %fontsize = fonts.ylabelfontsize;
  fontsize = fonts.ylabelfontsize+4;
  if(plot_struct.do_rescale_fluxes)
     h = ylabel(y_label,'fontsize',fontsize,'fontweight','bold');
  end
  annotate_plot(Xplotter,params,plot_struct,fonts, using_resistive_term, .27);
end

function plot_accum(plot_struct, t, E3, divPs_term, Ps13_x_term, Ps23_y_term, ...
  inertial_term, mu3i, ...
  using_resistive_term, resistive_term, residual, fluxes, params, species)

  fonts = plot_struct.fonts;
  % calculate accumulation integral of ohm's law terms
  %
  initial_yaxis_flux = get_initial_yaxis_flux(Xplotter(),params.outputdir);
  original_yaxis_flux = initial_yaxis_flux;
  int_E3             = integrate(t,E3);
  int_divPi_term     = integrate(t,divPs_term);
  int_Ps13_x_term    = integrate(t,Ps13_x_term);
  int_Ps23_y_term    = integrate(t,Ps23_y_term);
  int_inertial_term  = integrate(t,inertial_term);
  int_resistive_term = integrate(t,resistive_term);
  int_residual       = integrate(t,residual);

  if(plot_struct.do_rescale_times)
    t = t*params.ion_gyrofreq;
  end
  %scale_factor = 100./initial_yaxis_flux;
  %[scale_factor, y_label] = set_flux_scales(xplotter,plot_struct,params);
  [flux_scale_factor,initial_yaxis_flux] = get_flux_scale_factor(Xplotter(),...
    plot_struct,params,original_yaxis_flux);
  if(plot_struct.do_rescale_fluxes)
    int_E3             = flux_scale_factor * int_E3             ;
    int_divPi_term     = flux_scale_factor * int_divPi_term     ;
    int_Ps13_x_term    = flux_scale_factor * int_Ps13_x_term    ;
    int_Ps23_y_term    = flux_scale_factor * int_Ps23_y_term    ;
    int_inertial_term  = flux_scale_factor * int_inertial_term  ;
    mu3i               = flux_scale_factor * mu3i               ;
    int_resistive_term = flux_scale_factor * int_resistive_term ;
    int_residual       = flux_scale_factor * int_residual       ;
    if(fluxes.show_central_island_flux)
      fluxes.left_lost_flux = flux_scale_factor * fluxes.left_lost_flux;
    end
    if(fluxes.show_other_island_flux)
      fluxes.other_island_flux = flux_scale_factor * fluxes.other_island_flux;
    end
  end

  % plot cumulation integrals
  %
  figure(1+plot_struct.figure_shift);
  clf;
  hold on;
  if(plot_struct.do_rescale_times)
    xlabel(get_time_label(Xplotter(),plot_struct.do_rescale_times),...
      'fontsize',fonts.xlabelfontsize);
  end
  p1=plot(t,-int_E3             ,'b-' ); set(p1,'linewidth',4);
  p1=plot(t,-int_divPi_term     ,'r-'); set(p1,'linewidth',2);
  p1=plot(t,-mu3i               ,'k'  ); set(p1,'linewidth',1);
  p1=plot(t,-int_inertial_term  ,'k'); set(p1,'linewidth',2);
  if(using_resistive_term)
  p1=plot(t,-int_resistive_term ,'g'  ); set(p1,'linewidth',2);
  end
  p1=plot(t,-int_residual       ,'m'  ); set(p1,'linewidth',2);
  if(fluxes.show_central_island_flux)
  p1=plot(t,fluxes.left_lost_flux,'g-'); set(p1,'linewidth',2);
  end
  if(fluxes.show_other_island_flux)
  p1=plot(t,fluxes.other_island_flux,'m-'); set(p1,'linewidth',2);
  end
  p1=plot(t,-int_Ps13_x_term); set(p1,'linewidth',2,'color', [.6 .2 0],'linestyle','-');
  p1=plot(t,-int_Ps23_y_term); set(p1,'linewidth',2,'color', [.2 .6 0],'linestyle','-');
  hold off;

  % add title and axis labels
  %
  t1=title([ 'accumulation integral of "Ohm''s law" terms at the X-point ' ]);
  set(t1,'fontsize',fonts.titlefontsize,'FontWeight','bold');
  set(gca,'fontsize',fonts.ticlabelfontsize,'yaxislocation','left');
  set(gca,'FontWeight','bold');
  xlim([t(1), plot_struct.final_time]);
  %ylim([-4 40].*params.domain_scaling*params.E_0);
  if(plot_struct.show_percent_flux)
    ylim(plot_struct.ylims);
    %ylim([-2 10]);
  else
    ylim([-.2 .8].*initial_yaxis_flux);
  end
  axis on; box on; grid on;
  %set(gca,'xtick',0:100:500);
  %set(gca,'ytick',0:0.01:0.15);
  %set(gca,'plotboxaspectratio',[2 1 1]);
  set(gca,'linewidth',2);
  %
  % position the axes
  %
  axes_left=0.08;
  axes_bottom=0.12;
  axes_width = .84;
  axes_height = .78;
  set(gca,'position',[axes_left axes_bottom axes_width axes_height]);
  %
  % label the axes
  %
  h = xlabel(get_time_label(Xplotter(),plot_struct.do_rescale_times),...
    'fontsize',fonts.xlabelfontsize);
  position = get(h,'position');
  position(2)=-26;
  set(h,'position',position);
  %ylabel('percent of initial flux through y-axis');
  h = ylabel(get_y_label(plot_struct,params,initial_yaxis_flux),'fontsize', ...
    fonts.ylabelfontsize);
  position = get(h,'position');
  position(1)=-3.4;
  set(h,'position',position);
  %
  % annotate the plot with parameters
  %
  annotate_plot(Xplotter,params,plot_struct,fonts, using_resistive_term, .27);

  % create legend of quantities plotted
  %
  create_simple_legend(params,species,using_resistive_term,axes_left,fonts,...
    fluxes.show_central_island_flux, fluxes.show_other_island_flux);

end

function y_label = get_y_label(plot_struct,params, initial_yaxis_flux)
  % rescale fluxes if requested
  if(plot_struct.use_gauss)
    y_label = 'E_{gauss}';
  else
    y_label = 'E_{SI}';
  end
  y_inv_units_str = '';
  if(plot_struct.do_rescale_fluxes)
    if(plot_struct.use_gauss)
      y_units_str = 'B_0 v_A/(c\Omega_i)';
      y_inv_units_str = '\Omega_i c/(B_0 v_A)';
    else
      y_units_str = 'B_0 v_A/\Omega_i';
      y_inv_units_str = '\Omega_i /(B_0 v_A)';
    end
  end
  if(plot_struct.show_percent_flux)
    y_label = ['% of initial y-axis flux, 100% = ' sprintf('%5.4g',initial_yaxis_flux) '*' y_units_str];
  else
    y_label = [y_label ' ' y_units_str];
  end
end

function create_simple_legend(params,species,using_resistive_term,axes_left,fonts, ...
  show_central_island_flux, show_other_island_flux);

  if(species == 't')
    momentum_label = ['(m_iu_i-m_eu_e)_3/2'];
  else
    momentum_label = ['m_' species '(u_' species ')_3'];
  end

  qminus = get_qminus(species);
  labels = {...
     ... %'reconnected f{}l{}ux', ...
     'electric term', ...
     'pressure term', ...
     ... %['(u_' species ')_3 m_' species '/e'],...
     momentum_label, ...
     ... % 'inertia term', ... %['$\frac{m_' species '}{e} (u_' species ')_3$'],...
     'inertial term', ...
    };
  if(using_resistive_term)
    labels = [labels { '-$\int_0^t \hbox{resistive term}\ \ $' }];
  end
  labels = [labels {'residual term' }];
  if(show_central_island_flux)
    labels = [labels {'f{}l{}ux of center island'}];
  end
  if(show_other_island_flux)
    labels = [labels {'f{}l{}ux of side islands'}];
  end
  labels = [labels {'P_x_z_,_x term'}];
  labels = [labels {'P_y_z_,_y term'}];
  %legend_location = [axes_left+.09 .59  .155 .24 ];
  legend_location = [axes_left+.09 .59  .135 .24 ];
  %legend_location = 'northwest';
  l1=legend(labels,...
    ... % latex is screwed up; deletes certain characters (e.g. P,y,z)
    ... %'Interpreter', 'latex',
    'location', legend_location, ...
    'fontsize', fonts.legendfontsize);
end

function create_legend(params,species,using_resistive_term,axes_left,fonts, ...
  show_central_island_flux, show_other_island_flux);

  qminus = get_qminus(species);
  labels = {...
     '-$\int_0^t E_3$', ...
     [qminus '$\int_0^t \frac{(\nabla. P_' species ')_3}{e n_' species '}$'],...
     [qminus '$\frac{m_' species '}{e} (u_' species ')_3$'],...
     [qminus '$\int_0^t \frac{m_' species '}{e}\partial_t (u_' species ')_3$']...
    };
  if(using_resistive_term)
    labels = [labels { '-$\int_0^t \hbox{resistive term}\ \ $' }];
  end
  labels = [labels {'-$\int_0^t \hbox{residual}\ \ $' }];
  if(show_central_island_flux)
    labels = [labels {'f{}l{}ux of center island'}];
  end
  if(show_other_island_flux)
    labels = [labels {'f{}l{}ux of side islands'}];
  end
  legend_location = [axes_left+.08 .55  .12 .28 ];
  %legend_location = 'northwest';
  l1=legend(labels,...
    'Interpreter', 'latex', ...
    'location', legend_location, ...
    'fontsize', fonts.legendfontsize);
  %set(l1,'Interpreter','latex');
  %text('Interpreter','latex',...
  % 'String','$$\int_0^x\!\int_y dF(u,v)$$',...
  %  'Position',[.5 .5],...
  %   'FontSize',16);
  %set(l1,'location','northwest');
  %set(l1,'fontsize',fonts.legendfontsize);
end

function plot_resistivity(plot_struct, t, residual_slowing_rate, ...
  mu3i, M3i, params, species);

  fonts = plot_struct.fonts;
  % plot effective resistivity
  %
  figure(3+plot_struct.figure_shift);
  clf;
  hold on;
  p1=plot(t,residual_slowing_rate ,'k-'); set(p1,'linewidth',2);
  p1=plot(t,-mu3i ,'k-.'); set(p1,'linewidth',2);
  p1=plot(t,-M3i ,'r-.'); set(p1,'linewidth',2);
  hold off;
  labels = { ...
    'residual slowing rate', ...
    ['-m_' species ' (u_' species ')_3'], ...
    ['-(\rho_' species ' u_' species ')_3'], ...
  };
  l1 = legend(labels,...
    'location', 'northwest', ...
    'fontsize',fonts.legendfontsize);
  set(gca,'fontsize',fonts.ticlabelfontsize);
  t1=title(['effective residual slowing rate, i.e., (residual term)/(m_' ...
            species ' u_' species ')_3']);
  set(t1,'fontsize',fonts.titlefontsize);
  %axis([t(1) plot_struct.final_time -.1 .3]);
  axis 'auto x';
  %axis 'auto y';
  % xlim([t(1), plot_struct.final_time]);
  ylim([-.1 1].*params.domain_scaling);
  % axis 'auto y';
  axis on; box on; grid off;
  xlabel(get_time_label(Xplotter(),plot_struct.do_rescale_times),'fontsize', ...
    fonts.xlabelfontsize);
  hold off
  %set(gca,'xtick',0:100:500);
  %set(gca,'ytick',0:0.01:0.15);
  set(gca,'plotboxaspectratio',[1.5 1 1]);
end

