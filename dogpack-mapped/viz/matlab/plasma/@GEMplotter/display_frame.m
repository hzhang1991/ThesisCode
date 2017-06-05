function display_frame(c, time, data, var_name, varargin)
  s=c.s;

  % load mri
  % h = imshow(D(:,:,:,1)); set(h,'erasemode','xor');
  % for i = 2:150
  % set(h,'cdata',D(:,:,:,i))
  % pause
  % end
  
  % parse keyword arguments
  show_vectors = 1;
  raise_fig = 1; % 0;
  convect_fieldlines = c.s.plot_params.convect_fieldlines;
  color_range_args={};
  plot_type_args={};
  display_idx = 1;
  %
  idx=1;
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      case 'fig'; idx=idx+1;
        fig=varargin{idx};
        raise_fig=1;
      case 'display_string'; idx=idx+1;
        display_string=varargin{idx};
      case 'display_idx'; idx=idx+1;
        display_idx=varargin{idx};
      case 'color_range'; idx=idx+1;
        color_range_args=varargin(idx-1:idx);
        color_range = color_range_args{2};
        if(all(size(color_range)==[1, 2]))
          color_range = color_range';
        end
      case 'plot_type'; idx=idx+1;
        plot_type_args = varargin(idx-1:idx);
      case 'convect_fieldlines'; idx=idx+1;
        convect_fieldlines = varargin{idx};
      case 'show_vectors'; idx=idx+1;
        show_vectors = varargin{idx};
      otherwise
        %error(['invalid keyword: ' varargin{idx}]);
        % unrecognized keyword arguments are presumed to have one argument
        % and are ignored
        idx=idx+1;
    end
    idx=idx+1;
  end
  %
  % choose correct cell component if appropriate
  eval(command_to_get_var_from_cell_if_there('display_string','display_idx'));
  eval(command_to_get_var_from_cell_if_there('fig','display_idx'));

  % default keyword argument values
  %
  if(~exist('var_name','var'))
    var_name = 'var_name';
  end
  varIdx = -1;
  if(isfield(s.varMap.varIdx_s,var_name))
    varIdx = s.varMap.varIdx_s.(var_name);
  end
  if(~exist('display_string','var'))
    if(varIdx > -1)
      display_string = s.plot_params.display_names{varIdx};
    else
      display_string = var_name;
    end
  end

  if(exist('fig','var')&&~iscell(fig))
    use_default_fig=0;
  end
  if(use_default_fig)
    fig = varIdx;
  end
  % if(isa(fig,'cell'))
  if(fig>0)
    if(ishandle(fig) && ~raise_fig)
      set(0,'CurrentFigure',fig);
      set(fig,'erasemode','xor');
      % set(fig,'cdata',D(:,:,:,i))
    else
      figure(fig);
    end
  end
  clf;
  hold on;
  %
  plot_params = s.plot_params;
  %
  num_components = size(data,3);
  %num_components = s.varMap.size_s.(var_name);
  plot_grid = s.plot_grid;
  if(num_components==3)
   % display the third component with color
    if(numel(color_range_args)==0 && ~exist('color_range','var'))
      color_range = plot_params.color_ranges(:,varIdx);
      color_range_arguments = {'color_range', color_range};
    else
      color_range_arguments = color_range_args;
    end
     plot_scalar(c, data(:,:,3), ...
       'xflip', plot_params.flip.p(varIdx), ...
       'yflip', plot_params.flip.p(varIdx), ...
       color_range_arguments{:}, ...
       plot_type_args{:});
     %plot_scalar_cart2(c,data(:,:,3),plot_grid.xl,plot_grid.yl, ...
     %  plot_params.flip.p(varIdx),plot_params.flip.p(varIdx),...
     %  color_range,...
     %  plot_type_args{:});
  end
  if(num_components==2 || num_components==3)
     % use field lines to plot if vector field is supposed to be divergence-free
     do_plot_fieldlines = any(plot_params.noDiv(varIdx)==[1 3]);
     if(do_plot_fieldlines)
       % do not pass s.params.plot_dx; pass
       plot_fieldlines(c, plot_grid.xc,plot_grid.yc,data(:,:,1),data(:,:,2),...
         s.plotdims.dx,s.plotdims.dy,convect_fieldlines);
     end
     do_plot_vector = any(plot_params.noDiv(varIdx)==[2 3]);
     if(do_plot_vector && show_vectors)
       plot_vector(plot_grid.xc,plot_grid.yc,data(:,:,1),data(:,:,2),...
         plot_params.vector_scales(varIdx), ...
         plot_params.flip.p(varIdx), s.params.enforced_symmetry);
       % set_ylabel(c,rms);
     end
  elseif(num_components==1)
  % display with color
     plot_scalar(c, data,'var_name', var_name, color_range_args{:}, plot_type_args{:});
     %plot_scalar_cart2(data(:,:,1),plot_grid.xl,plot_grid.yl,... %time,var_name,...
     %   plot_params.flip.x(varIdx),plot_params.flip.y(varIdx),plot_params.color_ranges(:,varIdx));
  elseif(num_components==6)
     plot_eigs(data,plot_grid.xl,plot_grid.yl, ...
        plot_params.flip.p(varIdx),plot_params.color_ranges(:,varIdx));
  else
     error([' invalid num_components: ', ...
        num2str(plot_struct.size_s(varIdx))]);
  end
  use_long_labels = 0; %;~ ~ fig;
  set_axes(c,use_long_labels);
  set_xlabel(c,use_long_labels);
  set_title(c,display_string,time,use_long_labels);
  hold off;
end

function set_axes(c,use_long_labels)
    params = c.s.params;
    scaling = 1.;
    if(c.s.do_rescale) % nondimensionalize space units
      scaling = 1./params.ion_skindepth;
    end
    xhigh = params.xhigh * scaling;
    yhigh = params.yhigh * scaling;
    xlow = -xhigh;
    ylow = -yhigh;

    axis on; box on; grid off;
    axis('equal');
    axis([xlow xhigh ylow yhigh]);
    %xtick = (-12:4:12)*scaling;
    %ytick = (-6:2:6)*scaling;
    %set(gca,'xtick',xtick);
    %set(gca,'ytick',ytick);
    %disp(['c.s.do_rescale=' c.s.do_rescale]);
    %if(c.s.do_rescale)
      axes_left=0.05;
      axes_bottom=0.05;
      axes_width = .78;
      axes_height = .95;
    %else
    %  axes_left=0.05;
    %  axes_bottom=0.08;
    %  axes_width = .80;
    %  axes_height = .82;
    %end
    set(gca,...
      'fontsize', 16, ... %34,...
      'fontweight','bold',...
      'linewidth',2);
    %if(use_long_labels)
      set(gca,'position',[axes_left axes_bottom axes_width axes_height]);
    %end
end

% function set_ylabel(c,rms)
%   ylabel_str = ['rms hor. vec. mag. = ' sprintf('%.3g ', rms)];
%   h = ylabel(ylabel_str);
%   set(h,'fontsize', 20, 'fontweight','bold');
% end

function set_xlabel(c,use_long_labels)
    params = c.s.params;
    if(c.s.do_rescale)
      if(use_long_labels)
        sheet_width = num2str(params.sheet_thickness/params.ion_skindepth);
        xlabel_str = [...
          'spatial unit = ion inertial length = \delta_i' ...
          ', sheet width = ' num2str(sheet_width) '*\delta_i ' ];
      else
        xlabel_str = 'spatial unit = ion inertial length ';
      end
      if(~strcmp(params.model_name,'none'))
        xlabel_str = params.model_string;
      end
      h = xlabel(xlabel_str);
      set(h,'fontsize', 18, ... %36,...
        'fontweight','bold');
    end
end

function set_title(c,display_string,time,use_long_labels);
    params = c.s.params;
    if(c.s.do_rescale)
      time_str = ['t = ' num2str(time*params.ion_gyrofreq) ' / \Omega_i'];
    else
      time_str = ['t = ' num2str(time)];
    end
    if(c.s.do_rescale)
      iso_period_str = ['isotropization period = ' ...
        num2str(params.ion_iso_period*params.ion_gyrofreq) ' / \Omega_i'];
    else
      iso_period_str = ['isotropization period = ' ...
        num2str(params.ion_iso_period)];
    end
    if(params.ion_iso_period<0)
      iso_period_str = 'no isotropization';
    end
    if(use_long_labels)
      mesh_model_str = [' (' num2str(params.mx) 'x' num2str(params.my) ' grid, ' ...
        iso_period_str ')']; 
    else
      mesh_model_str = '';
    end

    t1 = title([display_string ' at ' time_str mesh_model_str ' ']); 
    set(t1,'fontsize', 18, ... %36,...
      'fontweight','bold');
end


% changed this to concatenate into one big array and plot in one
% command to avoid patchy-looking result from each quadrant
% having its own color map (in case color_range is [0,0]).
%
% xflip: -1 means reverse
% yflip: -1 means reverse
%
function plot_scalar_cart2(c,qvals,xl,yl,xflip,yflip, color_range, varargin);

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

    [mx,my]=size(qvals);
    % if I create a nonlinear scale I need to redefine
    % the numbers on the color bar.
    %
    % if(color_range(2)>0)
    %   qvals=color_range(2)*tanh(qvals_in/color_range(2));
    %   %qvals=color_range(2)*sign(qvals_in).*sqrt(tanh(abs(qvals_in)/color_range(2)));
    % else
    %   qvals=qvals_in;
    % end

    maxval=max(max(qvals));
    minval=min(min(qvals));
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
        %color_range=[-maxabs,maxabs];
        color_range=[minval;maxval];
      else
        color_range=[0,maxval];
      end
    end

    neg_x_is_virtual = c.s.params.neg_x_is_virtual;
    neg_y_is_virtual = c.s.params.neg_y_is_virtual;
    enforcing_rotational_symmetry = c.s.params.enforcing_rotational_symmetry;
    %
    % we expect size(xl)==size(yl)==[mx+1,my+1];
    % create arrays for plot coordinates and values.
    % (need to add unused upper row and column to accomodate
    % idiosyncrasy of contourf)
    %
    x_size = mx+1;
    y_size = my+1;
    sx = 0; % sx = starting x offset for lower (possibly right) corner
    sy = 0; % sy = starting y offset for (possibly upper) left corner
    if(neg_x_is_virtual)
      x_size = 2*mx+1;
      sx = mx;
    end
    if(neg_y_is_virtual)
      y_size = 2*my+1;
      sy = my;
    end
    %qplot=zeros(size(qvals)*2+1);
    qplot=zeros(x_size,y_size);
    %
    % lower left coordinates
    lx=zeros(x_size,y_size);
    ly=zeros(x_size,y_size);

    % copy into lower left
    if(neg_x_is_virtual && neg_y_is_virtual)
      assert(enforcing_rotational_symmetry==0);
      disp('display_frame.m: neg x and y are virtual')
      lx(1:mx,1:my)=-xl(mx+1:-1:2,my+1:-1:2);
      ly(1:mx,1:my)=-yl(mx+1:-1:2,my+1:-1:2);
      % should decide convention for not flipping (1?),
      % but this works as long as the user is consistent
      if(yflip==xflip)
        qplot(1:mx,1:my)=qvals(end:-1:1,end:-1:1);
      else
        qplot(1:mx,1:my)=-qvals(end:-1:1,end:-1:1);
      end
    end
    % copy into lower right
    if(neg_y_is_virtual)
      assert(enforcing_rotational_symmetry==0);
      lx(sx+1:end,1:my)= xl(:,my+1:-1:2);
      ly(sx+1:end,1:my)=-yl(:,my+1:-1:2);
      if(yflip==-1)
        % disp('flipping the y coordinate');
        qplot(sx+1:end-1,1:my)=-qvals(:,end:-1:1);
      else
        qplot(sx+1:end-1,1:my)= qvals(:,end:-1:1);
      end
    end
    % copy into upper left
    if(neg_x_is_virtual)
      lx(1:mx,sy+1:end)=-xl(mx+1:-1:2,:);
      ly(1:mx,sy+1:end)= yl(mx+1:-1:2,:);
      if(enforcing_rotational_symmetry)
        flip = xflip*yflip;
        qplot(1:mx,sy+1:end-1)=flip*qvals(end:-1:1,end:-1:1);
      else
        flip = xflip;
        qplot(1:mx,sy+1:end-1)= flip*qvals(end:-1:1,:);
      end
    end
    % copy into upper right
    %
    lx(sx+1:end,sy+1:end) = xl;
    ly(sx+1:end,sy+1:end) = yl;
    qplot(sx+1:end-1,sy+1:end-1)=qvals;

    switch lower(plot_type)
      case 'gradient' % schlieren
        [Ax,Ay] = gradient(qplot,s.plot_dims.dx,s.plot_dims.dy);
        grad_q = sqrt(Ax.*Ax + Ay.*Ay);
        phandle = surf(xc',yc',pgrad');
        shading flat;
        colormap(flipud(gray(2048)).^10);
      case 'pcolor'
        pcolor(lx,ly,qplot);
        % interp shading shifts the picture;
        % would need to use different coordinates
        %shading interp;
        shading flat;
        colormap('jet');
        if(exist('color_range','var'))
          caxis(color_range);
        end
        %else
          %caxis([0,mag]);
        %end
        c1=colorbar;
        set(c1,'fontsize',30);
      case 'contour'
        contourf(lx,ly,qplot,10); % show 10 discrete colors
        shading flat;
      otherwise
        error(['not supported: plot_type=' plot_type]);
    end

    set(gca,'fontsize',30);
end

function command = command_to_get_var_from_cell_if_there(vname,idx)
  command = ...
  ['if(~exist(''' vname ''',''var''));' ...
   'use_default_' vname ' = 1; return; else ' ...
   'use_default_' vname ' = 0; ' ...
    ' if(iscell(' vname ') && ' ...
         'numel(' vname ') >= ' num2str(idx) ');' ...
               '' vname ' ' ...
             '= ' vname '{' num2str(idx) '};' ...
   ' else ' ...
   'use_default_' vname ' = 2; ' ...
   ' end;' ...
   'end;'];
  % command = ...
  % ['if(exist(''' vname ''',''var'') && ' ...
  %      'iscell(' vname ') && ' ...
  %       'numel(' vname ') >= ' num2str(idx) ');' ...
  %             '' vname ' ' ...
  %           '= ' vname '{' num2str(idx) '};' ...
  %  'end;'];
end
%%%
