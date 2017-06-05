% printing examples:
%   plot_recon(mydir, final_time);
% to export Figure No. 2 to file recon.eps,
% with 100 dpi resolution, and using the EPS color graphics format, type
% print -f2 -r100 -depsc recon
%function plot_recon(xplotter,outputdir_in, final_time, do_rescale)
function plot_recon(xplotter,outputdir_in, varargin)

  setarg_outputdir;

  idx=1;
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      case 'final_time'; idx=idx+1;
        final_time=varargin{idx};
        idx=idx+1;
      case 'do_rescale'; idx=idx+1;
        do_rescale=varargin{idx}; idx=idx+1;
      otherwise
        error(['invalid keyword: ' varargin{idx}]);
    end
  end

  [fluxes] = get_fluxes(Xplotter(),outputdir);
  if(~exist('final_time','var'))
    final_time = fluxes.t(end);
  end
  plot_struct.final_time = final_time;

  params = get_GEM_params(outputdir);

  if(~exist('do_rescale','var')) do_rescale=1; end
  if(strcmp(params.model_name,'mhd'))
    do_rescale=0;
  end
  % first bit means rescale space-time.
  plot_struct.do_rescale_times = bitand(do_rescale,1);
  plot_struct.do_rescale_lengths = bitand(do_rescale,1);
  plot_struct.do_rescale_fluxes = bitand(do_rescale,2);
  plot_struct.use_gauss = bitand(do_rescale,4);
  if(plot_struct.do_rescale_times)
    plot_struct.final_time = final_time*params.ion_gyrofreq;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% read in the data and calculate quantities
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%
  %%% plot data
  %%%%%%%%%%%%%%%%%%%%

  plot_fluxes(fluxes, params, plot_struct);

  %%%%%%%%%%%%%%%%%%%%
  %%% print plot
  %%%%%%%%%%%%%%%%%%%%

  if(0) print_plot(params); end;

end

function plot_fluxes(fluxes, params, plot_struct)

  %fonts.ticlabelfontsize = 30;
  %fonts.annotationfontsize = 34;
  %fonts.xlabelfontsize = 32;
  %fonts.ylabelfontsize = 32;
  %fonts.titlefontsize = 40;
  %fonts.legendfontsize = 27;

  fonts.ticlabelfontsize = 20;
  fonts.annotationfontsize = 24;
  fonts.xlabelfontsize = 22;
  fonts.ylabelfontsize = 16;
  fonts.titlefontsize = 25;
  fonts.legendfontsize = 18;

  %fonts.ticlabelfontsize = 14;
  %fonts.annotationfontsize = 14;
  %fonts.xlabelfontsize = 14;
  %fonts.ylabelfontsize = 14;
  %fonts.titlefontsize = 18;
  %fonts.legendfontsize = 14;

  % plot fluxes
  %
  figure(100)
  clf
  hold on;

  % rescale per request
  %
  t = fluxes.t;
  if(plot_struct.do_rescale_times)
    t = t*params.ion_gyrofreq;
    xlabel(get_time_label(Xplotter(),plot_struct.do_rescale_times),...
      'fontsize',fonts.xlabelfontsize,'fontweight','bold');
  end
  [scale_factor, y_label] = set_flux_scales(Xplotter(),plot_struct,params);
  %scale_factor = 100./initial_yaxis_flux;
  ylabel(y_label,'fontsize',fonts.ylabelfontsize,'fontweight','bold');

  if(~test_equal(scale_factor,1))
    % fluxes.t  = fluxes.t*params.ion_gyrofreq;
    fluxes.change_in_bottom_to_rght_flux ...
      = fluxes.change_in_bottom_to_rght_flux * scale_factor;
    fluxes.rght_flux ...
      = fluxes.rght_flux * scale_factor;
    fluxes.left_flux ...
      = fluxes.left_flux * scale_factor;
    fluxes.recon_flux ...
      = fluxes.recon_flux * scale_factor;
    fluxes.left_lost_flux ...
      = fluxes.left_lost_flux * scale_factor;
    if(isfield(fluxes,'intE3'))
      fluxes.intE3 ...
        = fluxes.intE3 * scale_factor;
    end
  end

  p1=plot(t,fluxes.rght_flux,'-.','Color',[0,.5,.5]);  set(p1,'linewidth',2);
  p1=plot(t,fluxes.left_flux,'r-');  set(p1,'linewidth',2);
  p1=plot(t,fluxes.recon_flux,'k-'); set(p1,'linewidth',3);
  p1=plot(t,fluxes.change_in_bottom_to_rght_flux,'--','Color',[.5 0 .5]); set(p1,'linewidth',5);
  if(isfield(fluxes,'intE3'))
    p1=plot(t,-fluxes.intE3,'b-.');     set(p1,'linewidth',3);
  end
  if(fluxes.show_central_island_flux)
    % size of central island
    p1=plot(t,fluxes.left_lost_flux,'g-'); set(p1,'linewidth',4);
  end
  if(fluxes.show_other_island_flux)
    p1=plot(t,fluxes.other_island_flux,'m-'); set(p1,'linewidth',3);
  end
  %flux_check_rght_flux = fluxes.left_flux + fluxes.recon_flux;
  %p1=plot(t,flux_check_rght_flux,'y:'); % should agree with rght_flux
  hold off;
  set(gca,'fontsize',fonts.ticlabelfontsize);
  set(gca,'fontweight','bold');
  t1=title([ 'Reconnecting flux versus time']);
  set(t1,'fontsize',fonts.titlefontsize,'fontweight','bold');
  axis on; box on; grid on;
  set(gca,'linewidth',2);
  xlim([t(1), plot_struct.final_time]);
  ylim([-.2 1.2].*fluxes.left_flux(1)); % orig
  %ylim([-.2 1.2].*-fluxes.intE3(end)); % alt
  %
  set(gca,'linewidth',2);
  axes_left=0.08;
  axes_bottom=0.12;
  axes_width = .84;
  axes_height = .78;
  set(gca,'position',[axes_left axes_bottom axes_width axes_height]);

  annotate_plot(Xplotter(),params,plot_struct,fonts, 0, .02);
  % %model_problem_str = params.outputdir;
  % model_problem_string = get_model_problem_string(Xplotter(),params,...
  %   plot_struct.do_rescale_lengths);
  % limits = axis();
  % xpos = limits(1)+.02*(limits(2)-limits(1));
  % ypos = limits(3)+.04*(limits(4)-limits(3));
  % text(xpos, ypos, model_problem_string,'fontsize', ...
  %   fonts.annotationfontsize,'fontweight','bold');
  % xpos = limits(1)+.02*(limits(2)-limits(1));
  % ypos = limits(3)+.96*(limits(4)-limits(3));
  % method_string = get_method_string(Xplotter(),params);
  % text(xpos, ypos, method_string,'fontsize', ...
  %   fonts.annotationfontsize,'fontweight','bold');

  %set(gca,'plotboxaspectratio',[1.5 1 1]);
  labels = {...
     'flux exiting right side', ...
     'flux across y-axis',...
     'decreased y-axis flux',...
     'increased x-axis flux',
     };
  if(isfield(fluxes,'intE3'))
    labels = [labels {'-(integral)_0^t E_3(0)'}];
  end
  if(fluxes.show_central_island_flux)
    labels = [labels {'flux of center island'}];
  end
  if(fluxes.show_other_island_flux)
    labels = [labels {'flux of side islands'}];
  end
  %legend_location = [legend_left legend_bottom legend_width legend_height]
  %legend_location = 'west';
  legend_location = [axes_left+.18 (.5-.1)  .12 .21 ];
  %xpos = legend_location(1);
  %ypos = legend_location(2)+legend_location(4)+.4;
  %text(xpos, ypos, 'magnetic flux through boundaries of first quadrant of domain','fontsize',fonts.annotationfontsize,'fontweight','bold');
  l1=legend(labels, 'location', legend_location);
  %set(l1,'location',legend_location);
  %set(l1,'fontsize',legendfontsize);

  hold off;
end

function print_plot(params)

  cfl_str = ['_cfl=' num2str(params.cflv(2))];
  iso_period_str='';
  if(strcmp(params.model_name,'p10'))
    iso_period = num2str(params.ion_iso_period);
    if(params.ion_iso_period < 0)
      iso_period = 'infty';
    end
    iso_period_str = ['_iso_period=' iso_period];
  end

  basename = [getenv('HOME') '/figures'];
  filename = [basename 'recon_' ...
    params.model_name '_' num2str(params.mx) 'x' num2str(params.my) ...
    cfl_str iso_period_str];
  % replace periods with underscores so latex doesn't get confused.
  filename = regexprep(filename, '\.', '_');
  disp(['printing to file: ' filename]);
  print('-dpng',filename);
end

