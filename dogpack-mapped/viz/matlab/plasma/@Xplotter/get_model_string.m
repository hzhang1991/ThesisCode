function model_string = get_model_string(xplotter, params, do_rescale_times)
  model_name_string = get_model_name_string(params);
  %iso_period_string = get_iso_period_string(params,do_rescale_times);
  %if(numel(iso_period_string) > 0)
  %  model_string = [model_name_string ', ' iso_period_string];
  %else
    model_string = model_name_string;
  %end
end

function iso_period_string = get_iso_period_string(params, do_rescale_times)
  iso_period_string = '';
  iso_period = params.ion_iso_period;
  if(do_rescale_times)
    iso_period = iso_period * params.ion_gyrofreq;
    iso_period_str = [num2str(iso_period) '/\Omega_i'];
  else
    iso_period_str = [num2str(iso_period)];
  end

  if(strcmp(params.model_name,'p10'))
    assert(params.ion_iso_period==params.elc_iso_period);
    if(params.ion_iso_period<0)
      iso_period_string = 'no isotropization';
    else
      iso_period_string = ['isotropization period = ' iso_period_str];
    end
  elseif(strcmp(params.model_name,'g10'))
    iso_period_string = ['ion isotropization period = ' iso_period_str];
  end
end

function model_string = get_model_name_string(params)
  if(strcmp(params.model_name,'p05'))
    model_string = ['5-moment pair plasma'];
  elseif(strcmp(params.model_name,'p10'))
    model_string = ['10-moment pair plasma'];
  elseif(strcmp(params.model_name,'g05'))
    model_string = ['5-moment plasma'];
  elseif(strcmp(params.model_name,'g10'))
    model_string = ['10-moment plasma '];
  elseif(strcmp(params.model_name,'i10e5'))
    model_string = ['10-moment ions, 5-moment elcs '];
  elseif(strcmp(params.model_name,'mhd'))
    model_string = 'MHD';
  else
    error(['unsupported model: ' params.model_name]);
  end
end

