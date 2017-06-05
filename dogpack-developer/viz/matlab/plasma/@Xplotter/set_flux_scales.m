
function [scale_factor, y_label] = set_flux_scales(xplotter,plot_struct,params)
  % rescale fluxes if requested
  scale_factor = 1.;
  %if(plot_struct.use_gauss)
  %  y_label = 'E_{gauss}';
  %else
  %  y_label = 'E_{SI}';
  %end
  if(plot_struct.do_rescale_fluxes)
    scale_factor = scale_factor*params.ion_gyrofreq/params.E_0;
  end
  v0_str = 'v_A';
  if(plot_struct.use_gauss)
    y_label = ['units: B_0' v0_str '/(c\Omega_i)' ];
  else
    y_label = ['units: B_0' v0_str '/\Omega_i' ];
  end
end

