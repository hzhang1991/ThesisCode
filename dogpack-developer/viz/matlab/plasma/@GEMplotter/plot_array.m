
function plot_array(c,qplot,color_range,varargin)
  % parse keyword arguments
  plot_type='pcolor';
  %
  idx=1;
  while (idx<=numel(varargin))
    switch lower(varargin{idx})
      case 'plot_type'; idx=idx+1;
        plot_type = varargin{idx};
      otherwise
        idx=idx+1;
    end
    idx=idx+1;
  end

  if(~exist('color_range','var'))
    color_range=[0;0];
  end
  maxval=max(max(qplot));
  minval=min(min(qplot));
  if(~all(color_range==[0;0]))
    if(maxval>color_range(2))
      disp(['maximum value of ' num2str(maxval) ...
            ' exceeds colorbar max of ' num2str(color_range(2))]);
    end
    if(minval<color_range(1))
      disp(['minimum value of ' num2str(minval) ...
            ' falls below colorbar min of ' num2str(color_range(1))]);
    end
  else
    if(minval<0)
      %maxabs = max(maxval,-minval);
      %color_range=[-maxabs;maxabs];
      color_range=[minval;maxval];
    else
      color_range=[0;maxval];
    end
  end
  do_rescale_colors = 0;
  %if(any(color_range ~= [0;0]))
  %  do_rescale_colors = 1;
  %end
  if(do_rescale_colors)
    assert(color_range(2)>color_range(1));
    if(color_range(1)==0)
      cutoff_in=0;
      low = -color_range(2);
    else
      low = color_range(1);
      cutoff_in=low;
    end
    hgh = color_range(2);
    stretch=1.5;
    % map values to color values
    qcolors=remap(qplot,low,hgh,stretch);
    % get rounded y tick values
    [rounded_tick_vals,rounded_tick_vals_str] = ...
      get_tick_values(low,hgh,stretch,cutoff_in);
    extreme_colors=get_extreme_colors(low,hgh,stretch,cutoff_in);
    % map y tick values to color values
    % find colors of rounded ticks
    rounded_tick_colors=remap(rounded_tick_vals,low,hgh,stretch);
    % convert rounded_tick_vals_str to cell array
    rounded_tick_vals_str_arr=regexp(rounded_tick_vals_str,'(\S+)','tokens');
  else
    qcolors=qplot;
    extreme_colors=color_range;
  end
  switch plot_type
    case 'gradient' % schlieren
      [Ax,Ay] = gradient(qcolors,c.s.plotdims.dx,c.s.plotdims.dy);
      grad_q = sqrt(Ax.*Ax + Ay.*Ay);
      lx = c.s.plot_grid.lx;
      ly = c.s.plot_grid.ly;
      pcolor(lx,ly,grad_q);
      shading flat;
      colormap(flipud(gray(2048)).^10);
      %colormap('gray');
    case 'contour'
      cx = c.s.plot_grid.lx(1:end-1,1:end-1) + 0.5*c.s.plotdims.dx;
      cy = c.s.plot_grid.ly(1:end-1,1:end-1) + 0.5*c.s.plotdims.dy;
      contour_args = {20};
      if(any(color_range ~= [0;0]))
        caxis(extreme_colors);
        contour_levels = color_range(1)+(0:20)*((color_range(2)-color_range(1))/20);
        contour_args={contour_levels};
      end
      contourf(cx,cy,qcolors(1:end-1,1:end-1),contour_args{:},'k');
    case 'pcolor'
      lx = c.s.plot_grid.lx;
      ly = c.s.plot_grid.ly;
      % 'interp' is smoother than 'flat' but leaves whitespace
      % on the edges (but we could use zero-order extrapolation)
      shading_style = 'interp_full';
      if(strcmp(shading_style,'interp'))
        cx = lx(1:end-1,1:end-1) + 0.5*c.s.plotdims.dx;
        cy = ly(1:end-1,1:end-1) + 0.5*c.s.plotdims.dy;
        pcolor(cx,cy,qcolors(1:end-1,1:end-1));
        shading('interp');
      elseif(strcmp(shading_style,'interp_full'))
        plotdims = c.s.plotdims;
        cx1d = plotdims.dx*((-plotdims.pos_mx:plotdims.pos_mx+1)-0.5);
        cy1d = plotdims.dy*((-plotdims.pos_my:plotdims.pos_my+1)-0.5);
        [cx,cy]=ndgrid(cx1d,cy1d);
        qvc = zeros(size(qcolors)+1);
        % populate center
        qvc(2:end-1,2:end-1) = qcolors(1:end-1,1:end-1);
        %
        % extrapolate edges
        %
        qvc(2:end-1,1) = qcolors(1:end-1,1);
        qvc(2:end-1,end) = qcolors(1:end-1,end-1);
        qvc(1,2:end-1) = qcolors(1,1:end-1);
        qvc(end,2:end-1) = qcolors(end-1,1:end-1);
        %
        % extrapolate corners
        %
        qvc(1,1)=qcolors(1,1);
        qvc(1,end)=qcolors(1,end-1);
        qvc(end,1)=qcolors(end-1,1);
        qvc(end,end)=qcolors(end-1,end-1);
        %
        %max(max(qvc))
        %min(min(qvc))
        pcolor(cx,cy,qvc);
        shading('interp');
      else
        pcolor(lx,ly,qcolors);
        shading('flat');
      end
      colormap('jet');
      if(any(color_range ~= [0;0]))
        caxis(extreme_colors);
      end
    otherwise
      error(['invalid plot_type: ' plot_type]);
  end
  c1=colorbar;
  set(c1,'fontsize',30);
  if(do_rescale_colors)
    % convert cell array of cells to cell array of strings
    rounded_tick_vals_str_arr = [rounded_tick_vals_str_arr{:}];
    %rounded_tick_vals_str_arr = convert_arr_of_cells_to_arr_of_strings( ...
    %  rounded_tick_vals_str_arr);
    if(cutoff_in==low)
      rounded_tick_vals_str_arr{1}='-inf';
    end
    rounded_tick_vals_str_arr{end}='inf';
    set(c1,'ytick',rounded_tick_colors,'yticklabel',rounded_tick_vals_str_arr);
  end
  % This would be used for data with a periodic domain
  %colormap(hsv);
end

function arr = convert_arr_of_cells_to_arr_of_strings(arr)
  for i=1:numel(arr)
    arr(i)=arr{i};
  end
end

function [rounded_tick_vals,rounded_tick_vals_str] ...
  = get_tick_values(a,b,stretch_in,cutoff_in)

  stretch = stretch_in*.9999999999;
  ave=(a+b)*0.5;
  hdiff=(b-a)*0.5;
  stretch_hdiff = stretch*hdiff;
  % stretch about ave
  hghest = ave+stretch_hdiff;
  cutoff = ave+(cutoff_in-ave)*stretch;
  difference=(hghest-cutoff);
  tick_colors=cutoff+(0:.1:1)*difference;
  % evenly distribute ticks over range of colors
  % find values of ticks
  tick_vals = inverse_remap(tick_colors,a,b,stretch);
  % round values of ticks
  [rounded_tick_vals,rounded_tick_vals_str]=round_vals(tick_vals);
end

function [rounded_tick_vals,rounded_tick_vals_str]=round_vals(tick_vals)
  %ticks(end)
  % round to two digits
  rounded_tick_vals_str = num2str(tick_vals,'%.2g\n');
  rounded_tick_vals = str2num(rounded_tick_vals_str);
end

function color_extremes=get_extreme_colors(a,b,stretch,cutoff_in)
  ave=(a+b)*0.5;
  hdiff=(b-a)*0.5;
  stretch_hdiff = stretch*hdiff;
  % stretch about ave
  hghest = ave+stretch_hdiff;
  cutoff = ave+(cutoff_in-ave)*stretch;
  color_extremes=[cutoff,hghest];
end

function qplot=remap(qplot,a,b,stretch)
  ave=(b+a)*0.5;
  hdiff=(b-a)*0.5;
  rescale=hdiff*stretch;
  qplot = ave+rescale.*tanh((qplot-ave)./rescale);
end

function qplot=inverse_remap(qplot,a,b,stretch)
  stretch = 1.5;
  ave=(b+a)*0.5;
  hdiff=(b-a)*0.5;
  rescale=hdiff*stretch;
  qplot = ave+rescale.*atanh((qplot-ave)./rescale);
end

