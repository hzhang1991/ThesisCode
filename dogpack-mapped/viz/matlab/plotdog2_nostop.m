%
% points_per_dir = points per direction (spatial dimension)
%
% point_type = 1:   uniform points on each element
%            = 2:   Gauss-Legendre points on each element
%
function plotdog2_nostop(points_per_dir,outputdir_in,point_type)

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
  
  if (nargin>2)
    if(ischar(point_type))
      point_type = str2num(point_type);
    end
  else
    point_type=1;
  end
  
  if point_type~=1 && point_type~=2
    point_type = 1;
  end  
  
  if point_type==2 && points_per_dir>5
    points_per_dir=5;
  end

format long e;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FIND OUT IF CARTESIAN OR UNSTRUCTURED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fids  = fopen([outputdir,'/qhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/qhelp.dat  not found.']);
end
ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
if (ndims~=2)
    error(['Incorrect dimension, ndims must be 2. ndims = ',num2str(ndims)]);
end
GridType = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
GridTmp = GridType(1:9);
if (GridTmp=='Cartesian')
  GridType='Cartesian   ';
elseif (GridTmp=='Unstructu')
  GridType='Unstructured';
  error(['Currently only GridType=Cartesian is supported. GridType = ',GridType]);
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
  xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  ylow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  yhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  fclose(fids);
  
  % Grid information
  kmax = get_kmax(meth1,ndims);
  mx_old = mx;
  my_old = my;
  dx_old = (xhigh-xlow)/mx_old;
  dy_old = (yhigh-ylow)/my_old;
  mx = mx*points_per_dir;
  my = my*points_per_dir;
    
  if (point_type==1)
    xl = linspace(xlow,xhigh,mx+1);
    yl = linspace(ylow,yhigh,my+1);
    [xl,yl]=meshgrid(xl,yl);
    xl = xl';
    yl = yl';
  
    dx = (xhigh-xlow)/mx;
    dy = (yhigh-ylow)/my;
    xc = linspace(xlow+dx/2,xhigh-dx/2,mx);
    yc = linspace(ylow+dy/2,yhigh-dy/2,my);
    [xc,yc]=meshgrid(xc,yc);
    xc = xc';
    yc = yc';    
    
    dxi=1/points_per_dir;
    s1d =(-1+dxi):(2*dxi):(1-dxi);    

    kk=0;
    s2d = zeros(points_per_dir^2,2);    
    for jj=1:points_per_dir
      for ii=1:points_per_dir
        kk=kk+1;
        s2d(kk,1) = s1d(ii);
        s2d(kk,2) = s1d(jj);
      end
    end
    
  else

    sq3 = sqrt(3);
    sq5 = sqrt(5);
    sq7 = sqrt(7);
    
    if (points_per_dir==1)
      s1d = 0.0;
    elseif (points_per_dir==2)
      s1d = [-1.0/sq3; 1.0/sq3];
    elseif (points_per_dir==3)
      s1d = [-sq3/sq5; 0.0e0; sq3/sq5];
    elseif (points_per_dir==4)
      s1d = [-sqrt(3.0+sqrt(4.8))/sq7; -sqrt(3.0-sqrt(4.8))/sq7; ...
              sqrt(3.0-sqrt(4.8))/sq7;  sqrt(3.0+sqrt(4.8))/sq7];
    elseif (points_per_dir==5)
      s1d = [-sqrt(5.0 + sqrt(40.0/7.0))/3.0; ...
             -sqrt(5.0 - sqrt(40.0/7.0))/3.0; ...
              0.0; ...
              sqrt(5.0 - sqrt(40.0/7.0))/3.0; ...
              sqrt(5.0 + sqrt(40.0/7.0))/3.0];
    end

    kk=0;
    s2d = zeros(points_per_dir^2,2);
    for jj=1:points_per_dir
      for ii=1:points_per_dir      
        kk=kk+1;
        s2d(kk,1) = s1d(ii);
        s2d(kk,2) = s1d(jj);
      end
    end
    
    xc = zeros(mx,my);
    yc = zeros(mx,my);
        
    xline = (xlow+dx_old/2):dx_old:(xhigh-dx_old/2);
    yline = (ylow+dy_old/2):dy_old:(yhigh-dy_old/2);

    kk=1;
    for i=1:mx_old                 
      xx(kk:(kk+points_per_dir-1)) = xline(i)+(dx_old/2)*s1d;
      kk=kk+points_per_dir;
    end
    kk=1;
    for j=1:my_old                 
      yy(kk:(kk+points_per_dir-1)) = yline(j)+(dy_old/2)*s1d;
      kk=kk+points_per_dir;
    end

    [xc,yc]=meshgrid(xx,yy);
    xc=xc';
    yc=yc';
    
    xl = [xlow,(xx(1:end-1)+xx(2:end))/2,xhigh];
    yl = [ylow,(yy(1:end-1)+yy(2:end))/2,yhigh];
    [xl,yl]=meshgrid(xl,yl);
    xl = xl';
    yl = yl';
    
  end
  
  % Sample Legendre polynomial on the midpoint of each sub-element
  LegVals=GetCart2Legendre(kmax, s2d);
  
  xeps = max(0.015*(xhigh-xlow),0.015*(yhigh-ylow));
  yeps = xeps;
  
  q=-1;

  m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
              ' ) ? ']);
  disp(' ')
  if isempty(m)
    m=1;
  end

  for n1=0:nplot
          
      %% Solution -- q
      [q_data,time] = read_state2_cart(datafmt, outputdir, n1, 'q', ...
                                       mx_old, my_old, meqn, kmax, ...
                                       1:meqn);
      qsoln = sample_state2_cart(q_data, meth1, kmax, LegVals);
      clear q_data;
      qaug = zeros(mx+1,my+1,meqn);
      qaug(1:mx,1:my,1:meqn) = qsoln;
      
      if (maux>0)
        %% Aux variables -- aux
        [a_data,time] = read_state2_cart(datafmt, outputdir, n1, 'a', ...
                                         mx_old, my_old, maux, kmax, 1:maux);
        aux = sample_state2_cart(a_data, meth1, 1, LegVals);
        clear a_data;
        aux_aug = zeros(mx+1,my+1,maux);
        aux_aug(1:mx,1:my,1:maux) = aux;
      end
      
      % USER SUPPLIED FUNCTION: Plotting function
      plotq2_cart;
  end

end
