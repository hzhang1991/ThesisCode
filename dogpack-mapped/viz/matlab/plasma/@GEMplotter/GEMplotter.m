% g = GEMplotter(outputdir)
% plotout(GEMplotter(),outputdir,sample_rate)
% get_PlotParams
% plot_B(GEMplotter(outputdir,sample_rate),4)
%
% I am using GEMplotter to package the routines that are used by
% plotout to plot GEM problem frames.
%
% This is a mechanism to allow different user-accessible
% routines to access the same subroutines without exposing the
% subroutines to everyone.
%
function c = GEMplotter(outputdir_in, sample_rate, varargin)

  setarg_outputdir;
  convect_fieldlines = 0;

  idx=1;
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      case 'convect_fieldlines'; idx=idx+1;
        convect_fieldlines = varargin{idx};
      case 'do_rescale'; idx=idx+1;
        do_rescale = varargin{idx};
      otherwise
        error(['invalid keyword: ' varargin{idx}]);
    end
    idx=idx+1;
  end

  % There are two ways to handle the rescaling:
  % (1) rescale the output.  This is a little more flexible
  %     if we wish to use multiple rescalings in our output,
  %     as Bessho and Bhattacharjee do, and we only need
  %     to rescale what we choose to plot.
  % (2) rescale the input.  This seems the better way, since
  %     the input is limited and well-defined, whereas we might
  %     want to specify a plethora of outputs.
  % Of course both can be used on top of one another.
  % For now I just use (1).
  params = get_GEM_params(outputdir);

  if(~exist('sample_rate','var'))
    sample_rate = params.space_order;
  end
  if(~exist('do_rescale','var'))
    do_rescale=3;
  end
  if(strcmp(params.model_name,'mhd'))
    do_rescale=0;
  end

  % get plotting information
  %
  plotdims = get_plotdims(params,sample_rate,do_rescale);
  plot_grid = get_plot_grid(params, plotdims);
  
  % make structs of maps between components/variables and indices
  %
  tensorIdxMap = get_tensorIdxMap();
  stateVarMap = get_stateVarMap(params.model_name);
  auxVarMap = get_auxVarMap(params);
  varMap = get_varMap(params.model_name);
  assert_stateVarMap_agree_with_varMap(stateVarMap, varMap);
  assert_stateVarMap_agree_with_varMap(auxVarMap, varMap);

  % construct plot_params
  %
  plot_params = struct;
  %
  display_names = make_display_names(varMap, params);
  output_scales = make_output_scales(varMap, params);
  [color_ranges, vector_scales] = make_plotting_ranges(varMap, params);
  flip = make_flip(varMap);
  noDiv = make_noDiv(varMap);
  %
  plot_params.display_names=display_names;
  plot_params.output_scales = output_scales;
  plot_params.color_ranges   = color_ranges;
  plot_params.vector_scales  = vector_scales;
  plot_params.flip           = flip;
  plot_params.noDiv          = noDiv;
  plot_params.convect_fieldlines = convect_fieldlines;
  %
  stateInfo.state = [];
  stateInfo.time = -1;
  stateInfo.frame_number = -1;
  stateInfo.auxState = [];
  %
  s = struct;
  s.params = params;
  s.stateVarMap = stateVarMap;
  s.auxVarMap = auxVarMap;
  s.tensorIdxMap = tensorIdxMap;
  s.varMap = varMap;
  s.sample_rate = sample_rate;
  s.plot_params = plot_params;
  s.plotdims = plotdims;
  s.plot_grid = plot_grid;
  s.do_rescale = do_rescale;
  %
  s.stateInfo = stateInfo;
  %
  c.s = s;
  c = class(c,'GEMplotter');
end

%%% grid information

% get plot mesh parameters
function plotdims = get_plotdims(params,sample_rate,do_rescale)

  % sample_rate should be a positive or negative integer.
  assert(mod(sample_rate,1)==0);
  if(sample_rate > 0)
  % if sample_rate is postive then it is the number of samples per cell.
    % low/left edges of cells
    plotdims.dx = params.plot_dx/sample_rate;
    plotdims.dy = params.plot_dy/sample_rate;
    plotdims.mx = params.plot_mx*sample_rate;
    plotdims.my = params.plot_my*sample_rate;
  elseif(sample_rate < 0)
  % if sample_rate is negative then it is the number of cells
  % per (centered) sample.
    sample_period = -sample_rate;
    plotdims.dx = params.plot_dx*sample_period;
    plotdims.dy = params.plot_dy*sample_period;
    plotdims.mx = params.plot_mx/sample_period;
    plotdims.my = params.plot_my/sample_period;
  else
    error('invalid sample_rate')
  end
  % nondimensionalize space units
  if(do_rescale)
  plotdims.dx = plotdims.dx/params.ion_skindepth;
  plotdims.dy = plotdims.dy/params.ion_skindepth;
  plotdims.xlow = params.xlow/params.ion_skindepth;
  plotdims.xhigh = params.xhigh/params.ion_skindepth;
  plotdims.ylow = params.ylow/params.ion_skindepth;
  plotdims.yhigh = params.yhigh/params.ion_skindepth;
  end
  plotdims.sample_rate = sample_rate;

  % set derived parameters
  if(params.neg_x_is_virtual);
    plotdims.pos_mx = plotdims.mx;
  else
    plotdims.pos_mx = plotdims.mx/2;
  end; 
  if(params.neg_y_is_virtual);
    plotdims.pos_my = plotdims.my;
  else
    plotdims.pos_my = plotdims.my/2;
  end; 
  %plotdims.mx
  %plotdims.my
  %plotdims.pos_mx
  %plotdims.pos_my
end

function plot_grid=get_plot_grid(params, plotdims)
  [xl,yl,xl1d,yl1d]=get_left_grid2(params,plotdims);
  % centers of cells
  xc=xl(1:end-1,1:end-1)+(0.5*plotdims.dx);
  yc=yl(1:end-1,1:end-1)+(0.5*plotdims.dy);
  xc1d=xl1d(1:end-1)+(0.5*plotdims.dx);
  yc1d=yl1d(1:end-1)+(0.5*plotdims.dy);
  plot_grid.xl = xl;
  plot_grid.yl = yl;
  plot_grid.xc = xc;
  plot_grid.yc = yc;
  plot_grid.xc1d=xc1d;
  plot_grid.yc1d=yc1d;

  % coordinates for use by plot_scalar
  %
  % these are used by pcolor
  [lx,ly] = get_plot_grid2(params,xl,yl);
  % these are used by contourf
  %[cx,cy] = get_grid_centers2(xc,yc);
  plot_grid.lx = lx;
  plot_grid.ly = ly;
  %plot_grid.cx = cx;
  %plot_grid.cy = cy;
end

function [xl,yl,xl1d,yl1d]=get_left_grid2(params,plotdims)
  xl1d=plotdims.xlow+(0:plotdims.mx)*(plotdims.dx);
  %plotdims.ylow
  %plotdims.my*plotdims.dy
  yl1d=plotdims.ylow+(0:plotdims.my)*(plotdims.dy);
  % should the user do this?
  [xl,yl]=ndgrid(xl1d,yl1d);
end

function [cx,cy] = get_grid_centers2(xc,yc)
  mx = size(xc,1);
  my = size(xc,2);
  %[mx,my]=size(xl)
  cx=zeros(size(xc)*2);
  cy=zeros(size(yc)*2);
  cx(1:mx,1:my)=-xc(end:-1:1,end:-1:1);
  cy(1:mx,1:my)=-yc(end:-1:1,end:-1:1);
  cx(mx+1:end,1:my)= xc(:,end:-1:1);
  cy(mx+1:end,1:my)=-yc(:,end:-1:1);
  cx(1:mx,my+1:end)=-xc(end:-1:1,:);
  cy(1:mx,my+1:end)= yc(end:-1:1,:);
  cx(mx+1:end,my+1:end)= xc;
  cy(mx+1:end,my+1:end)= yc;
end

function [lx,ly] = get_plot_grid2(params,xl,yl)

  assert(all(size(xl)==size(yl)));

  mx = size(xl,1)-1;
  my = size(xl,2)-1;
  neg_x_is_virtual = params.neg_x_is_virtual;
  neg_y_is_virtual = params.neg_y_is_virtual;

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
  %
  % lower left coordinates
  lx=zeros(x_size,y_size);
  ly=zeros(x_size,y_size);

  % copy into lower left
  if(neg_x_is_virtual && neg_y_is_virtual)
    lx(1:mx,1:my)=-xl(mx+1:-1:2,my+1:-1:2);
    ly(1:mx,1:my)=-yl(mx+1:-1:2,my+1:-1:2);
  end
  % copy into lower right
  if(neg_y_is_virtual)
    lx(sx+1:end,1:my)= xl(:,my+1:-1:2);
    ly(sx+1:end,1:my)=-yl(:,my+1:-1:2);
  end
  % copy into upper left
  if(neg_x_is_virtual)
    lx(1:mx,sy+1:end)=-xl(mx+1:-1:2,:);
    ly(1:mx,sy+1:end)= yl(mx+1:-1:2,:);
  end
  % copy into upper right
  lx(sx+1:end,sy+1:end) = xl;
  ly(sx+1:end,sy+1:end) = yl;
end

%%% component information

function tensorIdxMap = get_tensorIdxMap();
  % map from pair of indices to linear index
  % (maybe not the best way to order indices of a symmetric matrix,
  % but entrenched in our code)
  tensorIdxMap.linearIdx = ...
    [[1,2,3]; ...
     [2,4,5]; ...
     [3,5,6]];
  % map from linear index to first index
  tensorIdxMap.linToIJ = ...
    [[1, 1];...
     [1, 2];...
     [1, 3];...
     [2, 2];...
     [2, 3];...
     [3, 3]];
  % example: return the 2nd index for the 5th linear index:
  % tensorIdxMap.linToIJ(5,2)
end

function c = get_stateVarMap(model_name)
  [name_arr, size_arr]=get_stateVar_names_and_sizes(model_name);
  c = make_stateVarMap(name_arr, size_arr);
end

function c = get_auxVarMap(params)
  [name_arr, size_arr]=get_auxVar_names_and_sizes(params);
  c = make_stateVarMap(name_arr, size_arr);
end

function assert_stateVarMap_agree_with_varMap(stateVarMap, varMap)
  for i=1:length(stateVarMap.name_arr)
    name = stateVarMap.name_arr{i};
    % This confirms that every component is a variable with the
    % same number of component indices
    assert(varMap.size_s.(name) == stateVarMap.size_s.(name));
  end
end

% This information should really be in the parameters.ini file, e.g. for p10:
%
%   varsize{6} = 1
%   varsize{5} = 1
%   varsize{4} = 2
%   varsize{3} = 6
%   varsize{2} = 3
%   varsize{1} = 1
%   var{6} = Psi
%   var{5} = E  
%   var{4} = B  
%   var{3} = Ni 
%   var{2} = Mi 
%   var{1} = rho_i
%
% or
%  ...
%  q{5} = nrg_i
%  q{4} = Mi_3
%  q{3} = Mi_2
%  q{2} = Mi_1
%  q{1} = rho_i
%   
% construct index mappings for components of state variable
%
function [name_arr, size_arr]=get_stateVar_names_and_sizes(model_name)
  if(strcmp(model_name,'mhd'))
    names_and_sizes = { ...
       'rho', 1, ...
       'M',    3, ...
       'energy', 1, ... % gas-dyamic energy plus magnetic energy
       'B',     3, ...
       'E',     3, ...
       'Psi',   1, ...
    };
  elseif(strcmp(model_name,'g10'))
    names_and_sizes = { ...
       'rho_i', 1, ...
       'Mi',    3, ...
       'Ni',    6, ...
       'rho_e', 1, ...
       'Me',    3, ...
       'Ne',    6, ...
       'B',     3, ...
       'E',     3, ...
       'Psi',   1, ...
       'phi',   1  ...
    };
  elseif(strcmp(model_name,'i10e5'))
    names_and_sizes = { ...
       'rho_i', 1, ...
       'Mi',    3, ...
       'Ni',    6, ...
       'rho_e', 1, ...
       'Me',    3, ...
       'nrg_e', 1, ...
       'B',     3, ...
       'E',     3, ...
       'Psi',   1, ...
       'phi',   1  ...
    };
  elseif(strcmp(model_name,'p10'))
    names_and_sizes = { ...
       'rho_i', 1, ...
       'Mi',    3, ...
       'Ni',    6, ...
       'B',     2, ...
       'E',     1, ...
       'Psi',   1  ...
    };
  elseif(strcmp(model_name,'g05'))
    names_and_sizes = { ...
       'rho_i', 1, ...
       'Mi',    3, ...
       'nrg_i', 1, ...
       'rho_e', 1, ...
       'Me',    3, ...
       'nrg_e', 1, ...
       'B',     3, ...
       'E',     3, ...
       'Psi',   1, ...
       'phi',   1, ...
    };
  elseif(strcmp(model_name,'p05'))
    names_and_sizes = { ...
       'rho_i', 1, ...
       'Mi',    3, ...
       'nrg_i', 1, ...
       'B',     2, ...
       'E',     1, ...
       'Psi',   1 ...
       'x_i',   1 ...
       'y_i',   1 ...
    };
  else
    error(['invalid model: ' model_name]);
  end
  name_arr = names_and_sizes(1:2:end);
  size_arr = cell2mat(names_and_sizes(2:2:end));
end

function [name_arr, size_arr]=get_auxVar_names_and_sizes(params)
  model_name = params.model_name;
  maux = params.maux;
  if(maux==0)
    name_arr = {};
    size_arr = {};
    return;
  end
  if(strcmp(model_name,'mhd'))
    names_and_sizes = { ...
       'aux1', 1, ...
    };
  elseif(strcmp(model_name,'g10'))
    names_and_sizes = { ...
       'aux1', 1, ...
       'aux2', 1  ...
    };
  elseif(strcmp(model_name,'i10e5'))
    names_and_sizes = { ...
       'aux1', 1, ...
       'aux2', 1  ...
    };
  elseif(strcmp(model_name,'p10'))
    names_and_sizes = { ...
       'aux1', 1, ...
    };
  elseif(strcmp(model_name,'g05'))
    names_and_sizes = { ...
       'aux1', 1, ...
       'aux2', 1  ...
    };
  elseif(strcmp(model_name,'p05'))
    names_and_sizes = { ...
       'aux1', 1, ...
    };
  else
    error(['invalid model: ' model_name]);
  end
  name_arr = names_and_sizes(1:2:end);
  size_arr = cell2mat(names_and_sizes(2:2:end));
end

% construct index mappings for all variables to be displayed
%
function varMap = get_varMap(model_name)
  Vn = 3;
  Bn = 3;
  sn = 1;
  if(is_symmetric(model_name))
    Vn = 1;
    Bn = 2;
    sn = 0;
  end
  var_names_and_sizes = { ...
     ... % MHD state variables
     'rho', 1, 'M', Vn, 'energy', 1, ...
     ... % other one-fluid gas state variables
     'N',  6,  'nrg',  1, ...
     ... % gas state variables
     'rho_i', 1, 'Mi', 3, 'Ni', 6, 'nrg_i', 1, ...
     'rho_e', 1, 'Me', 3, 'Ne', 6, 'nrg_e', 1, ...
     ... % electromagnetic state variables
     'B',Bn, 'E',Vn, 'Psi', 1, 'phi', 1, ...
     ... % tracking variables
     'x_i',1, 'y_i',1, ...
     'xi',1, 'yi',1, ...
     ... % primitive variables
     'ui', 3, 'Pi', 6, 'pi', 1, ...
     'ue', 3, 'Pe', 6, 'pe', 1, ...
     'u',  Vn, 'P', 6, 'p',  1, ...
     ... % spatial derivatives
     'divB', 1, ...
     'divE', 1, ...
     'divEc', 1, ...
     ... % charge densities
     'sgm_i', 1, ...
     'sgm_e', 1, ...
     'sgm',  sn, ...
     ... % currents
     'Ji', 3, ...
     'Je', 3, ...
     'J', Vn, ...
     'curlB', Vn, ...
     'JxB', Bn, ...
     'hall', 3, ...
     ... % Ohm's law terms
     'Bxu',Vn, 'Hal',Vn, 'pek',Vn, 'J_t',Vn, 'DFJ',Vn, 'dtJ',Vn, 'Eck',Vn, ...
     'nonideal',Vn, ...
     'resistivity',1, ...
     'Bxu3', 1, ...
     ... % and proxy Ohm's law terms
     'Bxui', 3, ...
     'Bxue', 3, ...
     'Bxui3', 1, ...
     'Bxue3', 1, ...
     'udotDu3_i', 1, ...
     'udotDu3_e', 1, ...
     'dtui3', 1, ...
     'dtue3', 1, ...
     'dtui', Vn, ...
     ... % components
     'ni', 1, ...
     'ne', 1, ...
     'Mi1', 1, ...
     'Mi2', 1, ...
     'Mi3', 1, ...
     'Me1', 1, ...
     'Me2', 1, ...
     'Me3', 1, ...
     'ui1', 1, ...
     'ui2', 1, ...
     'ui3', 1, ...
     'ue1', 1, ...
     'ue2', 1, ...
     'ue3', 1, ...
     'Ni11', 1, ...
     'Ni22', 1, ...
     'Ni33', 1, ...
     'Ni12', 1, ...
     'Ni13', 1, ...
     'Ni23', 1, ...
     'Ne11', 1, ...
     'Ne22', 1, ...
     'Ne33', 1, ...
     'Ne12', 1, ...
     'Ne13', 1, ...
     'Ne23', 1, ...
     'Pi11', 1, ...
     'Pi22', 1, ...
     'Pi33', 1, ...
     'Pi12', 1, ...
     'Pi13', 1, ...
     'Pi23', 1, ...
     'Pe11', 1, ...
     'Pe22', 1, ...
     'Pe33', 1, ...
     'Pe12', 1, ...
     'Pe13', 1, ...
     'Pe23', 1, ...
     'B1', 1, ...
     'B2', 1, ...
     'B3', 1, ...
     'E1', 1, ...
     'E2', 1, ...
     'E3', 1, ...
     ... % quantities derived from velocities 
     'ui1x', 1, ...
     'ue1x', 1, ...
     'ui2y', 1, ...
     'ue2y', 1, ...
     'divue', 1, ...
     'divui', 1, ...
     ... % quantities derived from pressure 
     'divPi3', 1, ...
     'divPe3', 1, ...
     'Pi13x', 1, ...
     'Pi23y', 1, ...
     'Pe13x', 1, ...
     'Pe23y', 1, ...
     'p', 1, ...
     'pi', 1, ...
     'pe', 1, ...
     'p1i', 1, ...
     'p2i', 1, ...
     'p3i', 1, ...
     'p1e', 1, ...
     'p2e', 1, ...
     'p3e', 1, ...
     'pbi', 1, ...
     'pbe', 1, ...
     'pmaxi', 1, ...
     'pmaxe', 1, ...
     'pmini', 1, ...
     'pmine', 1, ...
     'pratio_i', 1, ...
     'pratio_e', 1, ...
     'detPi', 1, ...
     'detPe', 1, ...
     'detTi', 1, ...
     'detTe', 1, ...
     'Ti', 1, ...
     'Te', 1, ...
     'T1i', 1, ...
     'T2i', 1, ...
     'T3i', 1, ...
     'T1e', 1, ...
     'T2e', 1, ...
     'T3e', 1, ...
     'entropy_i', 1, ...
     'entropy_e', 1, ...
     'Entropy_i', 1, ...
     'Entropy_e', 1, ...
     ... % wave speeds
     'maxspd_i', 1, ...
     'maxspd_e', 1, ...
     'cfi', 1, ...
     'cfe', 1, ...
     ... % braginskii closure
     'tau_i', 1, ...
     'tau_e', 1, ...
     'ptau_i', 1, ...
     'ptau_e', 1, ...
     'detPtau_i', 1, ...
     'detPtau_e', 1, ...
     'K_i', 1, ... % heat conductivity
     'K_e', 1, ... % heat conductivity
     'mu_i', 1, ... % viscosity
     'mu_e', 1, ... % viscosity
     'isorate_i', 1, ...
     'isorate_e', 1, ...
     'Sit', 1, ... %rate of production of entropy per volume for ions
     'Set', 1, ... %rate of production of entropy per volume for elcs
     'sit', 1, ... %rate of production of entropy per mass for ions
     'set', 1, ... %rate of production of entropy per mass for elcs
     'JdotEprime', 1, ...
     'aux1', 1, ...
     'aux2', 1, ...
     ... % time derivatives
     ... %'Mi_t', 'Me_t' ...
  };
  var_names = var_names_and_sizes(1:2:end);
  var_sizes = cell2mat(var_names_and_sizes(2:2:end));
  varMap = make_varMap(var_names, var_sizes);
end

function bool = is_symmetric(model_name)
  bool = strcmp(model_name,'p10')||strcmp(model_name,'p05');
end

% e.g.
%
% names_and_sizes = {...
%  'rho',  1, ...
%  'M',    3, ...
%  'N',    6};
%
% name_arr = names_and_sizes(1:2:end);
% size_arr = cell2mat(names_and_sizes(2:2:end));
% varMap = make_varMap(name_arr, size_arr);
%
function varMap = make_varMap(name_arr, size_arr);

  assert(iscell(name_arr));
  if(numel(size_arr)==0) size_arr = []; end
  assert(isnumeric(size_arr));
  assert(length(name_arr)==length(size_arr));

  size_s  = struct;
  varIdx_s = struct;

  for i=1:length(name_arr)
    the_name = name_arr{i};
    the_size = size_arr(i);
    size_s.(the_name) = the_size;
    varIdx_s.(the_name) = i;
  end

  % things to map between:
  %   name, idnum,
  % things to map to:
  %   size

  varMap = struct;
  varMap.name_arr = name_arr;
  varMap.size_arr = size_arr;
  varMap.varIdx_s = varIdx_s;
  varMap.size_s = size_s;
  varMap.num_variables = length(name_arr);
end

% e.g.
%
% names_and_sizes = {...
%  'rho',  1, ...
%  'M',    3, ...
%  'N',    6};
%
% name_arr = names_and_sizes(1:2:end);
% size_arr = cell2mat(names_and_sizes(2:2:end));
% stateVarMap = make_stateVarMap(name_arr, size_arr);
%
function stateVarMap = make_stateVarMap(name_arr, size_arr);

  stateVarMap = make_varMap(name_arr,size_arr);

  % add component index information to stateVarMap.

  lastStateIdx = 0;
  startStateIdx_s = struct;
  stateIndices_s = struct;
  for i=1:length(name_arr)
    the_name = name_arr{i};
    the_size = size_arr(i);
    start_idx = lastStateIdx+1;
    startStateIdx_s.(the_name) = start_idx;
    lastStateIdx = lastStateIdx + the_size;
    stateIndices_s.(the_name) = start_idx:lastStateIdx;
  end

  % Things to map among:
  %   name, idx,
  % Things to map to:
  %   
  %stateVarMap = struct; % we are adding to the struct, so we do NOT clear it.
  stateVarMap.startStateIdx_s = startStateIdx_s; % startStateIdx_s
  stateVarMap.stateIndices_s  = stateIndices_s;  % stateIndices_s
  stateVarMap.lastStateIdx    = lastStateIdx;    % lastStateIdx
end

% params is an argument only because I have in mind
% to put these names in the .ini file.
%
function display_names = make_display_names(varMap, params)
  % default display name is the variable name
  display_names = varMap.name_arr;
  varIdx_s = varMap.varIdx_s;

  scale = params.scale;
  display_names{varIdx_s.rho   } = 'total density';
  display_names{varIdx_s.M     } = 'total momentum';
  display_names{varIdx_s.energy} = 'total energy';
  %display_names{varIdx_s.rho_i } = 'ion density';
  display_names{varIdx_s.Mi    } = 'ion momentum';
  display_names{varIdx_s.Ni    } = 'ion energy tensor';
  display_names{varIdx_s.nrg_i } = 'ion scalar energy';
  %display_names{varIdx_s.rho_e } = 'electron density';
  display_names{varIdx_s.Me    } = 'electron momentum';
  display_names{varIdx_s.Ne    } = 'electron energy tensor';
  display_names{varIdx_s.nrg_e } = 'electron scalar energy';
  display_names{varIdx_s.rho   } = 'mass density';
  display_names{varIdx_s.M     } = 'net gas momentum';
  display_names{varIdx_s.N     } = 'net gas energy tensor';
  display_names{varIdx_s.nrg   } = 'net gas scalar energy';
  display_names{varIdx_s.B     } = 'magnetic field';
  display_names{varIdx_s.E     } = 'electric field';
  display_names{varIdx_s.Psi   } = 'magnetic correction potential';
  display_names{varIdx_s.phi   } = 'electric correction potential';
  display_names{varIdx_s.ui    } = 'ion fluid velocity';
  display_names{varIdx_s.Pi    } = 'ion pressure tensor';
  display_names{varIdx_s.pi    } = 'ion scalar pressure';
  display_names{varIdx_s.ue    } = 'electron fluid velocity';
  display_names{varIdx_s.Pe    } = 'electron pressure tensor';
  display_names{varIdx_s.pe    } = 'electron scalar pressure';
  display_names{varIdx_s.u     } = 'net fluid velocity';
  display_names{varIdx_s.P     } = 'net pressure tensor';
  display_names{varIdx_s.p     } = 'net scalar pressure';
  display_names{varIdx_s.divB  } = 'div(B)';
  display_names{varIdx_s.divE  } = 'div(E)';
  display_names{varIdx_s.divEc } = 'div(E)-\sigma/\epsilon_0';
  display_names{varIdx_s.sgm_i } = 'ion charge density';
  display_names{varIdx_s.sgm_e } = 'electron charge density';
  display_names{varIdx_s.sgm   } = 'net charge density';
  display_names{varIdx_s.Ji    } = 'ion current';
  display_names{varIdx_s.Je    } = 'electron current';
  display_names{varIdx_s.J     } = 'net current';
  display_names{varIdx_s.Bxu   } = 'B cross u';
  display_names{varIdx_s.Bxue  } = 'B cross u_e';
  display_names{varIdx_s.Bxui  } = 'B cross u_i';
  display_names{varIdx_s.Bxui3 } = '(B cross u_i)_3';
  display_names{varIdx_s.Bxue3 } = '(B cross u_e)_3';
  display_names{varIdx_s.dtui3 } = '(d_tu_i)_3';
  display_names{varIdx_s.dtue3 } = '(d_tu_e)_3';
  display_names{varIdx_s.Hal   } = 'Hall term';
  display_names{varIdx_s.pek   } = 'pressure term';
  %display_names{varIdx_s.J_t   } = 
  %display_names{varIdx_s.DFJ   } = 
  %display_names{varIdx_s.dtJ   } = 'inertial term';
  display_names{varIdx_s.dtui  } = 'inertial term';
  display_names{varIdx_s.nonideal  } = 'non-ideal electric field';
  display_names{varIdx_s.resistivity  } = 'effective resistivity';
  %display_names{varIdx_s.Eck   } = 
  display_names{varIdx_s.ni   } = 'n_i'; %'ion number density';
  display_names{varIdx_s.ne   } = 'n_e'; %'elc number density';
  %display_names{varIdx_s.Mi1   } = 
  %display_names{varIdx_s.Mi2   } = 
  %display_names{varIdx_s.Mi3   } = 
  %display_names{varIdx_s.Me1   } = 
  %display_names{varIdx_s.Me2   } = 
  %display_names{varIdx_s.Me3   } = 
  %display_names{varIdx_s.ui1   } = 
  %display_names{varIdx_s.ui2   } = 
  %display_names{varIdx_s.ui3   } = 
  %display_names{varIdx_s.ue1   } = 
  %display_names{varIdx_s.ue2   } = 
  %display_names{varIdx_s.ue3   } = 
  %display_names{varIdx_s.Ni11  } = 
  %display_names{varIdx_s.Ni22  } = 
  %display_names{varIdx_s.Ni33  } = 
  %display_names{varIdx_s.Ni12  } = 
  %display_names{varIdx_s.Ni13  } = 
  %display_names{varIdx_s.Ni23  } = 
  %display_names{varIdx_s.Ne11  } = 
  %display_names{varIdx_s.Ne22  } = 
  %display_names{varIdx_s.Ne33  } = 
  %display_names{varIdx_s.Pe12  } = 
  %display_names{varIdx_s.Pe13  } = 
  %display_names{varIdx_s.Pi11  } = 
  %display_names{varIdx_s.Pi22  } = 
  %display_names{varIdx_s.Pi33  } = 
  %display_names{varIdx_s.Pi12  } = 
  %display_names{varIdx_s.Pi13  } = 
  %display_names{varIdx_s.Pi23  } = 
  %display_names{varIdx_s.Pe11  } = 
  %display_names{varIdx_s.Pe22  } = 
  %display_names{varIdx_s.Pe33  } = 
  %display_names{varIdx_s.Pe12  } = 
  %display_names{varIdx_s.Pe13  } = 
  %display_names{varIdx_s.Pe23  } = 
  %display_names{varIdx_s.B1  } = 
  %display_names{varIdx_s.B2  } = 
  %display_names{varIdx_s.B3  } = 
  %display_names{varIdx_s.E1  } = 
  %display_names{varIdx_s.E2  } = 
  display_names{varIdx_s.E3  }    = 'E_3';
  display_names{varIdx_s.ui1x } = 'u_i_1_,_x';
  display_names{varIdx_s.ui2y } = 'u_i_2_,_y';
  display_names{varIdx_s.ue1x } = 'u_e_1_,_x';
  display_names{varIdx_s.ue2y } = 'u_e_2_,_y';
  display_names{varIdx_s.divue } = 'div(u_e)';
  display_names{varIdx_s.divui } = 'div(u_i)';
  display_names{varIdx_s.divPe3 } = 'div(P_e)_3';
  display_names{varIdx_s.divPi3 } = 'div(P_i)_3';
  %display_names{varIdx_s.Pi13x } = 
  %display_names{varIdx_s.Pi23y } = 
  %display_names{varIdx_s.Pe13x } = 
  %display_names{varIdx_s.Pe23y } = 
  %display_names{varIdx_s.p } = 
  %display_names{varIdx_s.pi } = 
  %display_names{varIdx_s.pe } = 
  display_names{varIdx_s.p1i } = 'p_1_i'; %'min ion pressure eigenvalue';
  display_names{varIdx_s.p2i } = 'p_2_i'; %'mid ion pressure eigenvalue';
  display_names{varIdx_s.p3i } = 'p_3_i'; %'max ion pressure eigenvalue';
  display_names{varIdx_s.p1e } = 'p_1_e'; %'min elc pressure eigenvalue';
  display_names{varIdx_s.p2e } = 'p_2_e'; %'mid elc pressure eigenvalue';
  display_names{varIdx_s.p3e } = 'p_3_e'; %'max elc pressure eigenvalue';
  display_names{varIdx_s.pbi } = 'b.P_i.b'; %'parallel ion pressure';
  display_names{varIdx_s.pbe } = 'b.P_e.b'; %'parallel elc pressure';
  %display_names{varIdx_s.pmaxi } = 
  %display_names{varIdx_s.pmaxe } = 
  %display_names{varIdx_s.pmini } = 
  %display_names{varIdx_s.pmine } = 
  %display_names{varIdx_s.pratio_i } = 
  %display_names{varIdx_s.pratio_e } = 
  display_names{varIdx_s.Ti } = 'T_i';
  display_names{varIdx_s.Te } = 'T_e';
  display_names{varIdx_s.K_i } = 'ion thermal conductivity';
  display_names{varIdx_s.K_e } = 'elc thermal conductivity';
  display_names{varIdx_s.mu_i } = 'ion viscosity';
  display_names{varIdx_s.mu_e } = 'elc viscosity';
  %display_names{varIdx_s.entropy_i } = 
  %display_names{varIdx_s.entropy_e } = 
  %display_names{varIdx_s.Entropy_i } = 
  %display_names{varIdx_s.Entropy_e } = 
  display_names{varIdx_s.udotDu3_i } = '(u_i.Grad(u_i))_3';
  display_names{varIdx_s.udotDu3_e } = '(u_e.Grad(u_e))_3';
  %display_names{varIdx_s.maxspd_i } = 
  %display_names{varIdx_s.maxspd_e } = 
  %display_names{varIdx_s.cfi   } = 
  %display_names{varIdx_s.cfe   } = 
  display_names{varIdx_s.tau_i } = 'Braginskii ion isotropization period';
  display_names{varIdx_s.tau_e } = 'Braginskii elc isotropization period';
  display_names{varIdx_s.detPtau_i } = 'detP ion isotropization period';
  display_names{varIdx_s.detPtau_e } = 'detP elc isotropization period';
  %display_names{varIdx_s.isoperiod_e } = 
  %display_names{varIdx_s.isorate_i } = 
  %display_names{varIdx_s.isorate_e } = 
  display_names{varIdx_s.JdotEprime } = 'MHD entropy production rate';
end

function output_scales = make_output_scales(varMap, params)
  output_scales = ones(1,varMap.num_variables);

  scale = params.scale;
  varIdx_s = varMap.varIdx_s;
  % This rescales all output by the Alfven speed.
  output_scales(:,varIdx_s.rho_i   ) = scale.rho_0;
  output_scales(:,varIdx_s.Mi      ) = scale.M_0;
  output_scales(:,varIdx_s.Ni      ) = scale.P_0;
  output_scales(:,varIdx_s.nrg_i   ) = scale.P_0;
  output_scales(:,varIdx_s.rho_e   ) = scale.rho_0;
  output_scales(:,varIdx_s.Me      ) = scale.M_0;
  output_scales(:,varIdx_s.Ne      ) = scale.P_0;
  output_scales(:,varIdx_s.nrg_e   ) = scale.P_0;
  output_scales(:,varIdx_s.rho     ) = scale.rho_0;
  output_scales(:,varIdx_s.M       ) = scale.M_0;
  output_scales(:,varIdx_s.N       ) = scale.P_0;
  output_scales(:,varIdx_s.nrg     ) = scale.P_0;
  output_scales(:,varIdx_s.B       ) = scale.B_0;
  output_scales(:,varIdx_s.E       ) = scale.E_0;
  output_scales(:,varIdx_s.Psi     ) = scale.E_0;
  output_scales(:,varIdx_s.phi     ) = scale.B_0;
  output_scales(:,varIdx_s.x_i     ) = scale.x_0*scale.rho_0;
  output_scales(:,varIdx_s.y_i     ) = scale.x_0*scale.rho_0;
  output_scales(:,varIdx_s.xi      ) = scale.x_0;
  output_scales(:,varIdx_s.yi      ) = scale.x_0;
  output_scales(:,varIdx_s.ui      ) = scale.u_0;
  output_scales(:,varIdx_s.Pi      ) = scale.P_0;
  output_scales(:,varIdx_s.pi      ) = scale.P_0;
  output_scales(:,varIdx_s.ue      ) = scale.u_0;
  output_scales(:,varIdx_s.Pe      ) = scale.P_0;
  output_scales(:,varIdx_s.pe      ) = scale.P_0;
  output_scales(:,varIdx_s.u       ) = scale.u_0;
  output_scales(:,varIdx_s.P       ) = scale.P_0;
  output_scales(:,varIdx_s.p       ) = scale.P_0;
  output_scales(:,varIdx_s.divB    ) = scale.B_0/scale.x_0;
  output_scales(:,varIdx_s.divE    ) = scale.E_0/scale.x_0;
  output_scales(:,varIdx_s.divEc   ) = scale.E_0/scale.x_0;
  output_scales(:,varIdx_s.sgm_i   ) = scale.n_0;
  output_scales(:,varIdx_s.sgm_e   ) = scale.n_0;
  output_scales(:,varIdx_s.sgm     ) = scale.n_0;
  output_scales(:,varIdx_s.Ji      ) = scale.J_0;
  output_scales(:,varIdx_s.Je      ) = scale.J_0;
  output_scales(:,varIdx_s.J       ) = scale.J_0;
  output_scales(:,varIdx_s.curlB   ) = scale.J_0;
  output_scales(:,varIdx_s.Bxu     ) = scale.E_0;
  output_scales(:,varIdx_s.Bxui    ) = scale.E_0;
  output_scales(:,varIdx_s.Bxue    ) = scale.E_0;
  output_scales(:,varIdx_s.Bxu3    ) = scale.E_0;
  output_scales(:,varIdx_s.Bxui3   ) = scale.E_0;
  output_scales(:,varIdx_s.Bxue3   ) = scale.E_0;
  output_scales(:,varIdx_s.JxB     ) = scale.J_0*scale.B_0;
  output_scales(:,varIdx_s.Hal     ) = scale.E_0;
  output_scales(:,varIdx_s.hall    ) = scale.E_0;
  output_scales(:,varIdx_s.pek     ) = scale.E_0;
  output_scales(:,varIdx_s.J_t     ) = scale.E_0;
  output_scales(:,varIdx_s.DFJ     ) = scale.E_0;
  output_scales(:,varIdx_s.dtJ     ) = scale.E_0;
  output_scales(:,varIdx_s.dtui3   ) = scale.E_0;
  output_scales(:,varIdx_s.dtue3   ) = scale.E_0;
  output_scales(:,varIdx_s.nonideal) = scale.E_0;
  output_scales(:,varIdx_s.resistivity) = scale.E_0/scale.J_0;
  output_scales(:,varIdx_s.Eck     ) = scale.E_0;
  output_scales(:,varIdx_s.ni      ) = scale.n_0;
  output_scales(:,varIdx_s.ne      ) = scale.n_0;
  output_scales(:,varIdx_s.Mi1     ) = scale.u_0;
  output_scales(:,varIdx_s.Mi2     ) = scale.u_0;
  output_scales(:,varIdx_s.Mi3     ) = scale.u_0;
  output_scales(:,varIdx_s.Me1     ) = scale.u_0;
  output_scales(:,varIdx_s.Me2     ) = scale.u_0;
  output_scales(:,varIdx_s.Me3     ) = scale.u_0;
  output_scales(:,varIdx_s.ui1     ) = scale.u_0;
  output_scales(:,varIdx_s.ui2     ) = scale.u_0;
  output_scales(:,varIdx_s.ui3     ) = scale.u_0;
  output_scales(:,varIdx_s.ue1     ) = scale.u_0;
  output_scales(:,varIdx_s.ue2     ) = scale.u_0;
  output_scales(:,varIdx_s.ue3     ) = scale.u_0;
  output_scales(:,varIdx_s.Ni11    ) = scale.P_0;
  output_scales(:,varIdx_s.Ni22    ) = scale.P_0;
  output_scales(:,varIdx_s.Ni33    ) = scale.P_0;
  output_scales(:,varIdx_s.Ni12    ) = scale.P_0;
  output_scales(:,varIdx_s.Ni13    ) = scale.P_0;
  output_scales(:,varIdx_s.Ni23    ) = scale.P_0;
  output_scales(:,varIdx_s.Ne11    ) = scale.P_0;
  output_scales(:,varIdx_s.Ne22    ) = scale.P_0;
  output_scales(:,varIdx_s.Ne33    ) = scale.P_0;
  output_scales(:,varIdx_s.Ne12    ) = scale.P_0;
  output_scales(:,varIdx_s.Ne13    ) = scale.P_0;
  output_scales(:,varIdx_s.Ne23    ) = scale.P_0;
  output_scales(:,varIdx_s.Pi11    ) = scale.P_0;
  output_scales(:,varIdx_s.Pi22    ) = scale.P_0;
  output_scales(:,varIdx_s.Pi33    ) = scale.P_0;
  output_scales(:,varIdx_s.Pi12    ) = scale.P_0;
  output_scales(:,varIdx_s.Pi13    ) = scale.P_0;
  output_scales(:,varIdx_s.Pi23    ) = scale.P_0;
  output_scales(:,varIdx_s.Pe11    ) = scale.P_0;
  output_scales(:,varIdx_s.Pe22    ) = scale.P_0;
  output_scales(:,varIdx_s.Pe33    ) = scale.P_0;
  output_scales(:,varIdx_s.Pe12    ) = scale.P_0;
  output_scales(:,varIdx_s.Pe13    ) = scale.P_0;
  output_scales(:,varIdx_s.Pe23    ) = scale.P_0;
  output_scales(:,varIdx_s.B1      ) = scale.B_0;
  output_scales(:,varIdx_s.B2      ) = scale.B_0;
  output_scales(:,varIdx_s.B3      ) = scale.B_0;
  output_scales(:,varIdx_s.E1      ) = scale.E_0;
  output_scales(:,varIdx_s.E2      ) = scale.E_0;
  output_scales(:,varIdx_s.E3      ) = scale.E_0;
  output_scales(:,varIdx_s.ui1x    ) = scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.ue1x    ) = scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.ui2y    ) = scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.ue2y    ) = scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.divui   ) = scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.divue   ) = scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.divPi3  ) = scale.P_0/scale.x_0;
  output_scales(:,varIdx_s.divPe3  ) = scale.P_0/scale.x_0;
  output_scales(:,varIdx_s.Pi13x   ) = scale.P_0/scale.x_0;
  output_scales(:,varIdx_s.Pi23y   ) = scale.P_0/scale.x_0;
  output_scales(:,varIdx_s.Pe13x   ) = scale.P_0/scale.x_0;
  output_scales(:,varIdx_s.Pe23y   ) = scale.P_0/scale.x_0;
  output_scales(:,varIdx_s.p       ) = scale.P_0;
  output_scales(:,varIdx_s.pi      ) = scale.P_0;
  output_scales(:,varIdx_s.pe      ) = scale.P_0;
  output_scales(:,varIdx_s.p1i     ) = scale.P_0;
  output_scales(:,varIdx_s.p2i     ) = scale.P_0;
  output_scales(:,varIdx_s.p3i     ) = scale.P_0;
  output_scales(:,varIdx_s.p1e     ) = scale.P_0;
  output_scales(:,varIdx_s.p2e     ) = scale.P_0;
  output_scales(:,varIdx_s.p3e     ) = scale.P_0;
  output_scales(:,varIdx_s.pbi     ) = scale.P_0;
  output_scales(:,varIdx_s.pbe     ) = scale.P_0;
  output_scales(:,varIdx_s.pmaxi   ) = scale.P_0;
  output_scales(:,varIdx_s.pmaxe   ) = scale.P_0;
  output_scales(:,varIdx_s.pmini   ) = scale.P_0;
  output_scales(:,varIdx_s.pmine   ) = scale.P_0;
  output_scales(:,varIdx_s.pratio_i) = scale.P_0;
  output_scales(:,varIdx_s.pratio_e) = scale.P_0;
  output_scales(:,varIdx_s.detPi   ) = scale.P_0^3;
  output_scales(:,varIdx_s.detPe   ) = scale.P_0^3;
  output_scales(:,varIdx_s.detTi   ) = scale.T_0^3;
  output_scales(:,varIdx_s.detTe   ) = scale.T_0^3;
  % E_t + q_x = 0, q=-K*T_x says E_t + (K*T_x)_x = 0
  scale_K_0 = scale.P_0*(scale.x_0*scale.u_0)/scale.T_0;
  output_scales(:,varIdx_s.K_i     ) = scale_K_0;
  output_scales(:,varIdx_s.K_e     ) = scale_K_0;
  % (rho*u)_t = stress_x, stress = mu*u_x says (rho*u)_t = mu*u_xx
  scale_mu_0 = scale.M_0*scale.x_0; % viscosity scale
  output_scales(:,varIdx_s.mu_i    ) = scale_mu_0;
  output_scales(:,varIdx_s.mu_e    ) = scale_mu_0;
  output_scales(:,varIdx_s.Ti      ) = scale.T_0;
  output_scales(:,varIdx_s.Te      ) = scale.T_0;
  output_scales(:,varIdx_s.T1i     ) = scale.T_0;
  output_scales(:,varIdx_s.T2i     ) = scale.T_0;
  output_scales(:,varIdx_s.T3i     ) = scale.T_0;
  output_scales(:,varIdx_s.T1e     ) = scale.T_0;
  output_scales(:,varIdx_s.T2e     ) = scale.T_0;
  output_scales(:,varIdx_s.T3e     ) = scale.T_0;
  output_scales(:,varIdx_s.entropy_i) = 1;
  output_scales(:,varIdx_s.entropy_e) = 1;
  output_scales(:,varIdx_s.Entropy_i) = 1;
  output_scales(:,varIdx_s.Entropy_e) = 1;
  output_scales(:,varIdx_s.udotDu3_i) = scale.u_0*scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.udotDu3_e) = scale.u_0*scale.u_0/scale.x_0;
  output_scales(:,varIdx_s.maxspd_i) = scale.u_0;
  output_scales(:,varIdx_s.maxspd_e) = scale.u_0;
  output_scales(:,varIdx_s.cfi     ) = scale.u_0;
  output_scales(:,varIdx_s.cfe     ) = scale.u_0;
  output_scales(:,varIdx_s.tau_i) = scale.t_0;
  output_scales(:,varIdx_s.tau_e) = scale.t_0;
  output_scales(:,varIdx_s.ptau_i) = scale.t_0;
  output_scales(:,varIdx_s.ptau_e) = scale.t_0;
  output_scales(:,varIdx_s.detPtau_i) = scale.t_0;
  output_scales(:,varIdx_s.detPtau_e) = scale.t_0;
  output_scales(:,varIdx_s.isorate_i) = 1./scale.t_0;
  output_scales(:,varIdx_s.isorate_e) = 1./scale.t_0;
  output_scales(:,varIdx_s.JdotEprime) = scale.J_0*scale.E_0;
end

function [color_ranges, vector_scales] = make_plotting_ranges(varMap, params)

  varIdx_s = varMap.varIdx_s;
  color_ranges=zeros(2,varMap.num_variables);
  %
  Em = 0.2;
  %
  color_ranges(:,varIdx_s.rho_i   ) = [0,0];
  color_ranges(:,varIdx_s.Mi      ) = [-0,0]; %[-0.2,.2];
  color_ranges(:,varIdx_s.Ni      ) = [0,0];
  color_ranges(:,varIdx_s.nrg_i   ) = [0,0];
  %
  color_ranges(:,varIdx_s.rho_e   ) = [0,0];
  color_ranges(:,varIdx_s.Me      ) = [-0,0];
  color_ranges(:,varIdx_s.Ne      ) = [-0,0];
  color_ranges(:,varIdx_s.nrg_e   ) = [0,0];
  %
  color_ranges(:,varIdx_s.rho     ) = [0,0];
  color_ranges(:,varIdx_s.M       ) = [-0,0];
  color_ranges(:,varIdx_s.N       ) = [-0,0];
  color_ranges(:,varIdx_s.nrg     ) = [0,1];
  %
  color_ranges(:,varIdx_s.B       ) = [params.B_guide-Em, params.B_guide+Em]; %[-Em,Em]; %[0., 2.];
  color_ranges(:,varIdx_s.E       ) = [-Em,Em];
  color_ranges(:,varIdx_s.Psi     ) = [-.05,.05];
  color_ranges(:,varIdx_s.phi     ) = [-0,0];
  %
  color_ranges(:,varIdx_s.Bxu     ) = [-Em,Em];
  color_ranges(:,varIdx_s.Hal     ) = [-Em,Em];
  color_ranges(:,varIdx_s.pek     ) = [-Em,Em];
  color_ranges(:,varIdx_s.J_t     ) = [-Em,Em];
  color_ranges(:,varIdx_s.DFJ     ) = [-Em,Em];
  color_ranges(:,varIdx_s.dtJ     ) = [-Em,Em];
  color_ranges(:,varIdx_s.dtui    ) = [-Em,Em];
  color_ranges(:,varIdx_s.nonideal) = [-Em,Em];
  color_ranges(:,varIdx_s.resistivity) = [-.1,.4];
  color_ranges(:,varIdx_s.Eck     ) = [-Em,Em];
  %
  color_ranges(:,varIdx_s.sgm     ) = [-1,1]*2e-1;
  %color_ranges(:,varIdx_s.divB    ) = [-.03,.03];
  %color_ranges(:,varIdx_s.divE    ) = [-10,10];  % why is this so large?
  %
  color_ranges(:,varIdx_s.u       ) = [-0,0];
  color_ranges(:,varIdx_s.ui2     ) = [-.2,.2];
  color_ranges(:,varIdx_s.J       ) = [-0,0];
  color_ranges(:,varIdx_s.Ji      ) = [-0,0];
  color_ranges(:,varIdx_s.Je      ) = [-0,0];
  color_ranges(:,varIdx_s.p       ) = [ 0,0];
  color_ranges(:,varIdx_s.P       ) = [ 0,0];
  color_ranges(:,varIdx_s.Entropy_i) = [ -5,1];
  color_ranges(:,varIdx_s.Entropy_e) = [ -5,1];
  color_ranges(:,varIdx_s.entropy_i) = [ -2,2];
  color_ranges(:,varIdx_s.entropy_e) = [ -2,2];
  %color_ranges(:,varIdx_s.aux1    ) = [ -2,2];
  color_ranges(:,varIdx_s.aux2    ) = [ -2,2];

  %color_ranges(:,varIdx_s.Pi12    ) = [-.006,.006];
  %color_ranges(:,varIdx_s.Pe12    ) = [-.006,.006];

  vector_scales = zeros(1,varMap.num_variables);

  vector_scales(:,varIdx_s.rho_i   ) = 0;
  vector_scales(:,varIdx_s.Mi      ) = .3; %color_ranges(2,varIdx_s.Mi);
  vector_scales(:,varIdx_s.Ni      ) = 0;
  vector_scales(:,varIdx_s.rho_e   ) = 0;
  vector_scales(:,varIdx_s.Me      ) = .3; %color_ranges(2,varIdx_s.Me);
  vector_scales(:,varIdx_s.Ne      ) = 0;
  vector_scales(:,varIdx_s.B       ) = 1;
  vector_scales(:,varIdx_s.E       ) = 1; %color_ranges(2,varIdx_s.E);
  vector_scales(:,varIdx_s.Psi     ) = 0;
  vector_scales(:,varIdx_s.phi     ) = 0;
  vector_scales(:,varIdx_s.Bxu     ) = color_ranges(2,varIdx_s.Bxu);
  vector_scales(:,varIdx_s.Hal     ) = color_ranges(2,varIdx_s.Hal);
  vector_scales(:,varIdx_s.pek     ) = color_ranges(2,varIdx_s.pek);
  vector_scales(:,varIdx_s.J_t     ) = color_ranges(2,varIdx_s.J_t);
  vector_scales(:,varIdx_s.DFJ     ) = color_ranges(2,varIdx_s.DFJ);
  vector_scales(:,varIdx_s.dtJ     ) = color_ranges(2,varIdx_s.dtJ);
  vector_scales(:,varIdx_s.dtui    ) = color_ranges(2,varIdx_s.dtui);
  vector_scales(:,varIdx_s.Eck     ) = color_ranges(2,varIdx_s.Eck);
  vector_scales(:,varIdx_s.sgm     ) = 0;
  %vector_scales(:,varIdx_s.divB    ) = 0;
  %vector_scales(:,varIdx_s.divE    ) = 0;
  vector_scales(:,varIdx_s.M       ) = .3;
  vector_scales(:,varIdx_s.u       ) = .3; %color_ranges(2,varIdx_s.u);
  vector_scales(:,varIdx_s.ui      ) = .3;
  vector_scales(:,varIdx_s.ue      ) = .3;
  vector_scales(:,varIdx_s.J       ) = color_ranges(2,varIdx_s.J);
end

function flip = make_flip(varMap);
  varIdx_s = varMap.varIdx_s;
  % mark pseudotensor quantities to be negated under reflection
  flip=struct;
  flip.p=ones(1,varMap.num_variables);
  flip.p([...
    varIdx_s.B,...
    varIdx_s.B3,...
    varIdx_s.Psi,...
    varIdx_s.divB...
    ]) = -1;
  flip.x=flip.p;
  flip.y=flip.p;
  % set non-tensor symmetries (e.g. components of tensors)
  flip.x([varIdx_s.Mi1]) = -1;
  flip.y([varIdx_s.Mi2]) = -1;
  flip.x([varIdx_s.ui1]) = -1;
  flip.y([varIdx_s.ui2]) = -1;
  flip.x([varIdx_s.Pi12]) = -1;
  flip.y([varIdx_s.Pi12]) = -1;
  flip.x([varIdx_s.Pe12]) = -1;
  flip.y([varIdx_s.Pe12]) = -1;
  flip.x([varIdx_s.Pi13]) = -1;
  flip.y([varIdx_s.Pi23]) = -1;
  flip.x([varIdx_s.Pe13]) = -1;
  flip.y([varIdx_s.Pe23]) = -1;
  flip.y([varIdx_s.B1]) = -1;
  flip.x([varIdx_s.B2]) = -1;
  flip.x([varIdx_s.x_i]) = -1;
  flip.y([varIdx_s.y_i]) = -1;
  flip.x([varIdx_s.xi]) = -1;
  flip.y([varIdx_s.yi]) = -1;
end

function noDiv = make_noDiv(varMap);
  % 0 = 4 = do not display vector field
  % 1 = represent vector field with contours
  % 2 = represent vector field with vectors
  % 3 = represent vector field both ways
  noDiv=ones(1,varMap.num_variables)*2;
  noDiv(varMap.varIdx_s.B) = 1;
end

function bool = is_symmetric_pair_model(c)
  model_name = c.s.params.model_name;
  bool = strcmp(model_name, 'p10') || strcmp(model_name, 'p05');
end
