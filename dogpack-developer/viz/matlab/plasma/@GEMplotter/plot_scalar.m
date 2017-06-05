function plot_scalar(c, qvals_in, varargin);

  % set keyword defaults
  %
  color_range=[0,0];

  % parse keyword arguments
  plot_type_args={};
  %
  idx=1;
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      %case 'qvals'; idx=idx+1;
      %  qvals=varargin{idx}; idx=idx+1;
      case 'var_name'; idx=idx+1;
        var_name=varargin{idx};
      case 'xflip'; idx=idx+1;
        xflip=varargin{idx};
      case 'yflip'; idx=idx+1;
        yflip=varargin{idx};
      case 'color_range';
        idx=idx+1;
        color_range=varargin{idx};
      case 'plot_type'; idx=idx+1;
        plot_type_args = varargin(idx-1:idx);
      otherwise
        error(['invalid keyword: ' varargin{idx}]);
    end
    idx=idx+1;
  end
  if(all(size(color_range)==[1, 2]))
    color_range = color_range';
  end

  % finish setting keyword defaults
  %
  if(exist('var_name','var'))
    varIdx_s = c.s.varMap.varIdx_s;
    if(isfield(varIdx_s,var_name))
      varIdx = varIdx_s.(var_name);
      title_str = c.s.plot_params.display_names{varIdx_s.(var_name)};
      if(~exist('color_range','var'))
        color_range = c.s.plot_params.color_ranges(:,varIdx);
      end
      if(~exist('xflip','var')); xflip = c.s.plot_params.flip.x(varIdx); end;
      if(~exist('yflip','var')); yflip = c.s.plot_params.flip.y(varIdx); end;
    else
      title_str = var_name;
    end
  end
  if(~exist('xflip','var')); xflip = 1; end
  if(~exist('yflip','var')); yflip = 1; end

  qplot = get_plotarray(c,qvals_in,xflip,yflip,1);

  plot_array(c,qplot,color_range,plot_type_args{:});

  % set default title string (callers might override).
  if(exist('title_str','var'))
    t1 = title(title_str);
    set(t1,'fontsize', 18, 'fontweight','bold');
  end
  set(gca,'fontsize',30);
end

% flip is plus or minus one
%
function qplot = get_plotarray(c,qvals,xflip,yflip,extra_cell)
  [mx,my]=size(qvals);
  % if(color_range(2)>0)
  %   qvals=color_range(2)*tanh(qvals/color_range(2));
  % end

  neg_x_is_virtual = c.s.params.neg_x_is_virtual;
  neg_y_is_virtual = c.s.params.neg_y_is_virtual;
  enforcing_rotational_symmetry = c.s.params.enforcing_rotational_symmetry;

  x_size = mx+extra_cell;
  y_size = my+extra_cell;
  sx = 0; % sx = starting x offset for 
  sy = 0; % sy = starting y offset of y
  if(neg_x_is_virtual)
    x_size = 2*mx+extra_cell;
    sx = mx;
  end
  if(neg_y_is_virtual)
    y_size = 2*my+extra_cell;
    sy = my;
  end
  %qplot=zeros(size(qvals)*2+1);
  qplot=zeros(x_size,y_size);

  % copy into lower left
  if(neg_x_is_virtual && neg_y_is_virtual)
    assert(enforcing_rotational_symmetry==0);
    if(yflip==xflip)
      qplot(1:mx,1:my)=qvals(end:-1:1,end:-1:1);
    else
      qplot(1:mx,1:my)=-qvals(end:-1:1,end:-1:1);
    end
  end
  % copy into lower right
  if(neg_y_is_virtual)
    assert(enforcing_rotational_symmetry==0);
    if(yflip==-1)
      % disp('flipping the y coordinate');
      qplot(sx+1:sx+mx,1:my)=-qvals(:,end:-1:1);
    else
      qplot(sx+1:sx+mx,1:my)= qvals(:,end:-1:1);
    end
  end
  % copy into upper left
  if(neg_x_is_virtual)
    if(enforcing_rotational_symmetry)
      flip = xflip*yflip;
      qplot(1:mx,sy+1:sy+my)=flip*qvals(end:-1:1,end:-1:1);
    else
      flip = xflip;
      qplot(1:mx,sy+1:sy+my)= flip*qvals(end:-1:1,:);
    end
  end
  % copy into upper right
  qplot(sx+1:sx+mx,sy+1:sy+my)=qvals;
end

