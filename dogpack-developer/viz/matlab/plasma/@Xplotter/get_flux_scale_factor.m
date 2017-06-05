
function [flux_scale_factor,initial_yaxis_flux] = get_flux_scale_factor(xplotter,...
  plot_struct,params,initial_yaxis_flux)

  flux_scale_factor = 1.;
  if(plot_struct.do_rescale_times)
    flux_scale_factor = flux_scale_factor * params.ion_gyrofreq;
  end
  if(plot_struct.do_rescale_fluxes)
    flux_scale_factor = flux_scale_factor / params.E_0;
    initial_yaxis_flux = flux_scale_factor * initial_yaxis_flux ;

    if(plot_struct.show_percent_flux)
      flux_scale_factor = 100. * flux_scale_factor / initial_yaxis_flux;
    end
  end
end
