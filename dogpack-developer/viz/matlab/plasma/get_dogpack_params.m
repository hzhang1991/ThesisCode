function params = get_dogpack_params(outputdir, params)
  ini_fileName = [outputdir '/parameters.ini'];
  if(~exist('params')); params=struct; end;
  params = read_more_params_from_ini(ini_fileName, params, ...
    ... % primary method parameters
      'space_order', ... % space discretization
      'time_order', 'cflv',      ... % time discretization
    ... % other dogpack parameters
      'nout', 'tfinal', 'dtv', 'nv',  ...
      'meqn', 'maux', 'mbc', ...
      'datafmt' ...
    );
  out_ini_fileName = [outputdir '/out_parameters.ini'];
  fid = fopen(out_ini_fileName,'r');
  if(fid~=-1) fclose(fid); end
  if(fid==-1)
    % for backward compatibility with data
    params = read_more_params_from_ini(ini_fileName, params, ...
        'mx', 'my', 'xlow', 'xhigh', 'ylow', 'yhigh', 'dx','dy');
    params.out_ini = 0;
    params.plot_mx = params.mx;
    params.plot_my = params.my;
  else
    % new way
    params.out_ini = 1;
    params = read_more_params_from_ini(out_ini_fileName, params, ...
        'meqn', 'mx', 'my', ...
        'xlow', 'xhigh', 'ylow', 'yhigh', ...
        'dx', 'dy', ...
        'plot_mx','plot_my','plot_dx','plot_dy'...
      );
    % hack for compatibility with brief data window missing mx and my
    if(~isfield(params,'mx'))
      params = read_more_params_from_ini(ini_fileName, params, ...
        'mx', 'my');
    end
  end

end
