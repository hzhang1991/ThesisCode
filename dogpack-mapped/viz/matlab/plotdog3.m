function plotdog3(points_per_dir, outputdir_in)
%PLOTDOG3    Generic 3D plotting routine for DoGPack.
%
% PLOTDOG3(points_per_dir, outputdir_in) is the main plotting routine used in
% DoGPack for all of the 3D problems.
%
% Input parameters:
%
%   points_per_dir = points per grid element.  Default = 1.
%
%   outputdir_in - string identifying the output directory.  This can be
%   relative or an absolute pathname.  Default = 'output'.
%
%   point_type = 1:   uniform points on each element
%              = 2:   Gauss-Legendre points on each element
%
% Output parameters:
%
%    None.
%
% See also: plotq3_cart, plotdog1, plotdog2, plotdog4

  % Parse the input parameters (set default values for points_per_dir,
  % point_type and output folder name:
  global outputdir
  parse_input_args;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FIND OUT IF CARTESIAN OR UNSTRUCTURED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fids  = fopen([outputdir,'/qhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/qhelp.dat  not found.']);
end
ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
if (ndims~=3)
    error(['Incorrect dimension, ndims must be 3. ndims = ',num2str(ndims)]);
end
GridType = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
GridTmp = GridType(1:9);
if (GridTmp=='Cartesian')
  GridType='Cartesian   ';
elseif (GridTmp=='Unstructu')
  GridType='Unstructured';
  error(['In 3D currently only GridType=Cartesian is supported. GridType = ',GridType]);
else
  error(['Incorrect GridType, must be Cartesian or Unstructured. GridType = ',GridType]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp(['        GridType = ',GridType]);
disp(['  points_per_dir = ',num2str(points_per_dir)]);
disp(['      point_type = ',num2str(point_type)]);
disp(['       outputdir = ',outputdir]);

disp(' ');

if (GridType=='Cartesian   ')
    
  meqn    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  maux    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  nplot   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  meth1   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  datafmt = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  %
  mx      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  my      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  mz      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  ylow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  yhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  zlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  zhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  fclose(fids);
  
  % Grid information
  kmax = get_kmax(meth1,ndims);
  mx_old = mx;
  my_old = my;
  mz_old = mz;
  dx_old = (xhigh-xlow)/mx_old;
  dy_old = (yhigh-ylow)/my_old;
  dz_old = (zhigh-zlow)/mz_old;
  mx = mx*points_per_dir;
  my = my*points_per_dir;
  mz = mz*points_per_dir;
    
  xl = linspace(xlow,xhigh,mx+1);
  yl = linspace(ylow,yhigh,my+1);
  zl = linspace(zlow,zhigh,mz+1);
  [yl,xl,zl]=meshgrid(yl,xl,zl);
  
  dx = (xhigh-xlow)/mx;
  dy = (yhigh-ylow)/my;
  dz = (zhigh-zlow)/mz;
  xc = linspace(xlow+dx/2,xhigh-dx/2,mx);
  yc = linspace(ylow+dy/2,yhigh-dy/2,my);
  zc = linspace(zlow+dz/2,zhigh-dz/2,mz);
  [yc,xc,zc]=meshgrid(yc,xc,zc);
  
  dxi=1/points_per_dir;
  s1d =(-1+dxi):(2*dxi):(1-dxi);    
  
  mm=0;
  s3d = zeros(points_per_dir^3,2);    
  for kk=1:points_per_dir
      for jj=1:points_per_dir
          for ii=1:points_per_dir
              mm=mm+1;
              s3d(mm,1) = s1d(ii);
              s3d(mm,2) = s1d(jj);
              s3d(mm,3) = s1d(kk);
          end
      end
  end
    
  % Sample Legendre polynomials at s3d:
  LegVals = GetCart3Legendre(kmax, s3d);
  
  xeps = max([0.015*(xhigh-xlow),0.015*(yhigh-ylow),0.015*(xhigh-xlow),0.015*(zhigh-zlow)]);
  yeps = xeps;
  zeps = xeps;
  
  q=-1;

  m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
              ' ) ? ']);
  disp(' ')
  if isempty(m)
    m=1;
  end

  kn = 0;

  n = 0;
  nf = 0;
  n1 = -1;

  while(nf~=-1)
    nf  = input([ ' Plot which frame ( 0 - ',num2str(nplot),...
                  ' ) [type -1 or q to quit] ? ']);
    if isempty(nf)
      n1 = n1 + 1;
      nf = 0;
    else
      n1 = nf;
    end
    if n1> nplot
      disp(' ');
      disp(' End of plots ');
      disp(' ');
      n1 = nplot;
    end
    if (nf~=-1)
      
      %% Solution -- q
      [q_data,time] = read_state3_cart(datafmt, outputdir, n1, 'q', ...
                                       mx_old, my_old, mz_old, meqn, kmax, ...
                                       1:meqn);
      qsoln = sample_state3_cart(q_data, meth1, kmax, LegVals);
      clear q_data;
      qaug = zeros(mx+1,my+1,mz+1,meqn);
      qaug(1:mx,1:my,1:mz,1:meqn) = qsoln;
      
      if (maux>0)
        %% Aux variables -- aux
        [a_data,time] = read_state3_cart(datafmt, outputdir, n1, 'a', ...
                                         mx_old, my_old, mz_old, maux, ...
                                         kmax, 1:maux);

        % TODO - I can't find the cart_mod function in the vizualization
        % library (-DS)
        % aux = sample_state3_cart_mod(a_data, meth1, kmax, LegVals);
        aux = sample_state3_cart(a_data, meth1, kmax, LegVals);
        clear a_data;
        aux_aug = zeros(mx+1,my+1,mz+1,maux);
        aux_aug(1:mx,1:my,1:mz,1:maux) = aux;
      end
      
      % USER SUPPLIED FUNCTION: Plotting function
      plotq3_cart;
    end

  end
  disp(' ')

    
else
  disp(' ');
  disp([' Error in plotdog3.m: GridType = ',GridType,' is not ' ...
                      'supported.']);
  disp(' ');
end

