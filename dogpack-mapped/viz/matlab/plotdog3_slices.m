function plotdog3_slices(points_per_dir,outputdir_in)
%PLOTDOG3(points_per_dir, outputdir_in)
%
% points_per_dir = points per direction (spatial dimension)
%                  (default: points_per_dir = 1)
%
% outputdir      = directory to access for DoGPack output data
%                  (default: outputdir = 'output')
%
  global outputdir
  outputdir='output';
  
  if(nargin<1)
    points_per_dir = 1;
  end 
  
  if (ischar(points_per_dir))
    points_per_dir = str2num(points_per_dir);
  end

  if(nargin>1)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end

  if points_per_dir<1
    points_per_dir = 1;
  end
  
  point_type = 1;

format long e;


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
  
  % Read-in slice information from qhelp_slice.dat
  fids  = fopen([outputdir,'/qhelp_slice.dat'],'r');
  if fids==-1
      error(['File  ',outputdir,'/qhelp_slice.dat  not found.']);
  end
  numslices = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  if (numslices<1)
      error(['Number of slices must be at least 1, numslices = ',num2str(numslices)]);
  end
  slicekind = zeros(numslices,1);
  sliceidx  = zeros(numslices,1);
  for ns=1:numslices
      slicekind(ns) = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
      sliceidx(ns)  = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  end
  
  disp(['       numslices = ',num2str(numslices)]);
  for ns=1:numslices
        disp(['    slicekind(',num2str(ns),') = ',num2str(slicekind(ns)),...
              ', sliceidx(',num2str(ns),') = ',num2str(sliceidx(ns))]);
  end
  disp(' ');
  
  % Grid information
  kmax = get_kmax(meth1,2);
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
  
  xc_z(1:mx,1:my) = xc(:,:,1);
  yc_z(1:mx,1:my) = yc(:,:,1);
  
  xc_y(1:mx,1:mz) = xc(:,1,:);
  zc_y(1:mx,1:mz) = zc(:,1,:);

  yc_x(1:my,1:mz) = yc(1,:,:);
  zc_x(1:my,1:mz) = zc(1,:,:);
  
  xl_z(1:mx+1,1:my+1) = xl(:,:,1);
  yl_z(1:mx+1,1:my+1) = yl(:,:,1);
  
  xl_y(1:mx+1,1:mz+1) = xl(:,1,:);
  zl_y(1:mx+1,1:mz+1) = zl(:,1,:);
  
  yl_x(1:my+1,1:mz+1) = yl(1,:,:);
  zl_x(1:my+1,1:mz+1) = zl(1,:,:);
  
  dxi=1/points_per_dir;
  s1d =(-1+dxi):(2*dxi):(1-dxi);    
  
  mm=0;
  s2d = zeros(points_per_dir^2,2);    
  for jj=1:points_per_dir
      for ii=1:points_per_dir
          mm=mm+1;
          s2d(mm,1) = s1d(ii);
          s2d(mm,2) = s1d(jj);
      end
  end
    
  % Sample Legendre polynomial on the midpoint of each sub-element
  LegVals=GetCart2Legendre(kmax, s2d);
  
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
     
      % dimension each to be able to store each slice
      mmax = max([mx,my,mz]);
      qaug = zeros(mmax+1,mmax+1,meqn,numslices);
      if (maux>0)
          aux_aug = zeros(mmax+1,mmax+1,maux,numslices);
      end
      ms1 = zeros(numslices,1);
      ms2 = zeros(numslices,1);
      
      %% Loop over each slice
      for ns=1:numslices
          %% Get dimensions
          switch slicekind(ns)
            case 1
              m1 = mx;
              m2 = my;
              m1_old = mx_old;
              m2_old = my_old;
            case 2
              m1 = mx;
              m2 = mz;
              m1_old = mx_old;
              m2_old = mz_old;
            case 3
              m1 = my;
              m2 = mz;
              m1_old = my_old;
              m2_old = mz_old;
          end
          
          ms1(ns) = m1;
          ms2(ns) = m2;
          
          %% Solution -- q
          basefilename = [outputdir '/' 'q' ...
                          sprintf('%04d', n1) ...
                         '_slice' sprintf('%04d', ns)];
          [q_data,time] = read_state2_cart(datafmt, outputdir, n1, 'q', ...
                                           m1_old, m2_old, meqn, kmax, ...
                                           1:meqn, basefilename);
          qsoln = sample_state2_cart(q_data, meth1, kmax, LegVals);
          clear q_data;
          qaug(1:m1,1:m2,1:meqn,ns) = qsoln(1:m1,1:m2,1:meqn);
          clear qsoln;
      
          if (maux>0)
              %% Aux variables -- aux
              basefilename = [outputdir '/' 'a' ...
                              sprintf('%04d', n1) ...
                              '_slice' sprintf('%04d', ns)];
              [a_data,time] = read_state2_cart(datafmt, outputdir, n1, 'a', ...
                                               m1_old, m2_old, maux, ...
                                               kmax, 1:maux, basefilename);
              aux = sample_state2_cart(a_data, meth1, kmax, LegVals);
              clear a_data;
              aux_aug(1:m1,1:m2,1:maux,ns) = aux(1:m1,1:m2,1:maux);
              clear aux;
          end
      end
      
      % USER SUPPLIED FUNCTION: Plotting function
      plotq3_cart_slices;
    end

  end
  disp(' ')

    
else
  disp(' ');
  disp([' Error in plotdog3_slices.m: GridType = ',GridType,' is not ' ...
                      'supported.']);
  disp(' ');
end

