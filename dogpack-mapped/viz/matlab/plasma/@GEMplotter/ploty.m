
% plot a slice of a state variable along the y-axis
%
% usage: ploty(g,{'E','pek','dtui','Bxu'},34)
%        ploty(g,'E',34)
%
% example: show a movie of y profiles
%
% g=GEMplotter(outputdir,1)
% for i=0:80; ploty(g,{'-E','-pek','-dtui','-Bxu'},i,'ylim',[-.1,.2]); pause(.3); end
% for i=0:80;ploty(g,{'-E', '-pek', '-dtui', '-Bxu'},i,'ylim',[-.1,.2],'axis','x'); pause(.3); end
% ploty(g,{'-2*E','-.4*J','-2*nonideal','resistivity'},34,'ylim',[-.05,.4], 'var_name', 'resistivity');
%
% args = {'title', 'resistivity profile at time %time / \Omega_i', 'xlabel', 'y-axis (distance per skin depth)', 'ylabel', 'units nondimensionalized by B_0, n_0, and v_A'};
% ploty(g,{'resistivity','-.4*J','-2*nonideal','-2*E'},34,'ylim',[-.05,.45], args{:});
%
function ploty(c,eval_strings,frame_number,varargin)

  % set default values of keyword arguments
  do_print=0;
  %
  % parse keyword arguments
  %
  % set defaults
  axis_str = 'y';
  title_str='%axis_str-axis profile at time %time / \Omega_i';
  % could get this from var_name.
  %ylabel_str = { ...
  %  '-(electric field)/(B_0 v_A)', ...
  %  '-J_z/(e n_0 v_A)' };
  ylabel_str={'(nondimensionalized units)',''};
  xlabel_str = ''; %'distance per skin depth';
  component=1;
  %
  unused_args = cell(size(varargin));
  num_unused_args = 0;
  idx=1;
  two_scales=0;
  legend_labels={};
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      case 'title'; idx=idx+1;
        title_str=varargin{idx};
      case 'ylabel'; idx=idx+1;
        ylabel_str{1}=varargin{idx};
      case 'ylabel2'; idx=idx+1;
        ylabel_str{2}=varargin{idx};
      case 'xlabel'; idx=idx+1;
        xlabel_str=varargin{idx};
      case 'legend'; idx=idx+1; % legend labels
        legend_labels=varargin{idx};
      case 'ylim'; idx=idx+1;
        ylims=varargin{idx};
      case 'axis'; idx=idx+1;
        axis_str=varargin{idx};
      case 'component'; idx=idx+1;
        component=varargin{idx};
      case '-print';
        do_print=1; 
      case '-two';
        two_scales=1; 
      otherwise; idx=idx+1;
        % all unused arguments are assumed to have
        % a corresponding value
        num_unused_args = num_unused_args + 2;
        unused_args{num_unused_args-1} = varargin{idx-1};
        unused_args{num_unused_args} = varargin{idx};
    end
    idx=idx+1;
  end
  unused_args = unused_args(1:num_unused_args);
  %
  % check arguments
  assert(axis_str=='x' || axis_str=='y');
  % convert eval_strings to cell array if string
  if(ischar(eval_strings))
    eval_strings={eval_strings};
  end
  [state,time] = read_state(c,frame_number);
  % reduce the amount of data we have to process
  % (but plotting takes much more time than reading and calculation anyway)
  % assuming enforced symmetry
  if(axis_str=='y')
    state = state(1:3,:,:,:);
  else
    state = state(:,1:3,:,:);
  end
  c.s.stateInfo.state = state;
  c.s.stateInfo.time = time;
  %
  % get slice values
  %
  slice_vals=cell(size(eval_strings));
  labels=cell(size(eval_strings));
  for v=1:numel(eval_strings)
    % get this slice
    %
    eval_string=eval_strings{v};
    % this is only second-order accurate.
    % I should ideally generate the legendre coefficients
    % in the cells where I need them and then
    % evaluate them at the points where I need them,
    % but I don't want to rewrite all the sample functions
    % to output Legendre coefficients.  As an alternative
    % I could extrapolate.  But this is robust.
    [var_vals,var_name,display_string] ...
      = sample_expr_from_state(c,eval_string,unused_args{:});
    %[var_vals,var_name,display_string] ...
    %  = sample_expr_from_state(c,eval_string,state);
    if(numel(legend_labels) >= v)
      labels{v} = legend_labels{v};;
    else
      labels{v} = display_string;
    end
    if(axis_str=='y')
      slice = var_vals(1,:,component);
    else
      slice = var_vals(:,1,component)';
    end
    %
    % apply symmetries of var_name
    %
    varIdx = c.s.varMap.varIdx_s.(var_name);
    if(axis_str=='y')
      flip = c.s.plot_params.flip.y(varIdx);
    else
      flip = c.s.plot_params.flip.x(varIdx);
    end
    assert(flip*flip==1);
    reflected_slice = slice(end:-1:1);
    if(flip==-1) reflected_slice = -reflected_slice; end
    full_slice = [reflected_slice, slice];
    slice_vals{v} = full_slice;
  end
  %title_str = regexprep(title_str,'%time',num2str(time));
  %title_str = regexprep(title_str,'%axis_str',axis_str);
  repl_expressions = {'%axis_str', '%time'};
  repl_values = {axis_str,num2str(time)};
  title_str = regexprep(title_str, repl_expressions, repl_values);

  % move this to get_plot_grid in GEMplotter.m?
  if(axis_str=='y')
    pc1d = c.s.plot_grid.yc1d;
  else
    pc1d = c.s.plot_grid.xc1d;
  end
  c1d = [-pc1d(end:-1:1) pc1d];
  xlims=[c1d(1), c1d(end)];

  % plot the slices
  %markers={'b-','r-','k','g','m'};
  colors=    { 'b', 'r', 'k', 'g', 'm'};
  linestyles={ '-', '-', '-', '-', '-'};
  linewidths={ 5  ,   4,   3,   3,   3};
  clf;
  hold off;
  hold on;
  if(two_scales && numel(slice_vals)==2)
    %plotyy(c1d,slice_vals{1},c1d,slice_vals{2},'linewidth',linewidths{1:2});
    H=1:2;
    [AX,H(1),H(2)]=plotyy(c1d,slice_vals{1},c1d,slice_vals{2},'plot');
    for idx=1:2
      y_label_handle = get(AX(idx),'Ylabel');
      set(y_label_handle,...
        'String',ylabel_str{idx},...
         'color',colors{idx}, ...
        'fontsize',24);

      set(H(idx),...
          'color',colors{idx}, ...
          'linestyle',linestyles{idx}, ...
          'linewidth',linewidths{idx});
      set(AX(idx),'ycolor',colors{idx});
      if(exist('ylims','var'))
        set(AX(idx),'ylim',ylims{idx});
      end
      %set(AX(idx),'Xlim',xlims);
    end
    set(AX,...
       'fontsize',14,...
       'fontweight','bold',...
       'linewidth',2,...
       'box','on');
    set(AX,'Xlim',xlims);
    %
    % position = [left bottom width height]
    %position=get(AX(1),'position')
    %position(3) = position(3)*.88;
    %position(1) = (1-position(3))/2;
    %set(AX,'position',position);
    %
    %curr_ylim = get(AX,'Ylim')
    %curr_ylim{1} = curr_ylim{1}*1.5;
    %curr_ylim{2} = curr_ylim{2}*1.5/2;
    %set(AX(1),'ylim',curr_ylim{1});
    %set(AX(2),'ylim',curr_ylim{2});
  else
    for v=1:numel(slice_vals)
      vmod = mod(v,numel(colors));
      if(~vmod); vmod = numel(colors); end;
      h = plot(c1d,slice_vals{v});
      set(h, ...
        'color',colors{vmod}, ...
        'linestyle',linestyles{vmod}, ...
        'linewidth',linewidths{vmod});
    end
  legend(labels,'location','northwest','fontsize',12); %16
  title(title_str,'fontsize',26,'fontweight','bold');
  set(gca,'fontsize',20,'fontweight','bold','linewidth',2,'box','on');
  if(exist('ylims','var')) ylim(ylims); end;
  ylabel(ylabel_str{1},'fontsize',22,'fontweight','bold');
  xlim(xlims);
  end
  grid on;
  if(numel(xlabel_str)~=0)
    xlabel(xlabel_str);
  end
  hold off;
  if(do_print)
    basename = [getenv('HOME') '/figures'];
    filename = [basename '/' [eval_strings{:}] '_' num2str(frame_number) '.eps'];
    disp(['printing to file: ' filename]);
    print('-depsc',filename);
  end
end

%function minus = strip_minus(eval_string)
%    minus_label='';
%    minus=1;
%    if(eval_string(1)=='-')
%      minus=-1;
%    end
%end
