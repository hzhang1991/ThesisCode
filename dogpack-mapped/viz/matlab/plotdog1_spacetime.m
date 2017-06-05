% THIS FUNCTION WILL CURRENTLY NOT WORK.
%
% Right now, it doesn't set up the correct points to be evaluated.  The call
% to SampleBasis1 will fail.  See plotdog1.m
%
% I don't see any applications that call this function.  Can we remove it?
% (-DS).
%
function plotdog1_spacetime(points_per_dir,outputdir_in)

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

format long e;

disp(' ');
disp(['       outputdir = ',outputdir]);
disp(['  points_per_dir = ',num2str(points_per_dir)]);
disp(' ');

fids  = fopen([outputdir,'/qhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/qhelp.dat  not found.']);
end
nplot = fscanf(fids,'%d',1);
meqn  = fscanf(fids,'%d',1);
maux  = fscanf(fids,'%d',1);
meth1 = fscanf(fids,'%d',1);
mx    = fscanf(fids,'%d',1);
xlow  = fscanf(fids,'%e',1);
xhigh = fscanf(fids,'%e',1);
dx    = fscanf(fids,'%e',1);
fclose(fids);

% Grid information
mx_old = mx;
mx = mx*points_per_dir;
dx = (xhigh-xlow)/mx;
xc = transpose(linspace(xlow+dx/2,xhigh-dx/2,mx));

% Sample basis functions on mesh
% size of phi = (points_per_dir,meth1);
phi = SampleBasis1(points_per_dir,meth1);
  
q=-1;

m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
            ' ) ? ']);
disp(' ')
if isempty(m)
  m=1;
end

qtrack = zeros(mx,nplot+1,meqn);
auxtrack = zeros(mx,nplot+1,maux);
for n1=0:nplot
  %% Solution -- q
  fname = [outputdir,'/',num2str(n1+10000),'.dat'];
  fname(8) = 'q';	      
  fids = fopen(fname,'r');
  time(n1+1,1) = fscanf(fids,'%e',1);
  qtmp = fscanf(fids,'%e',[1,inf]);
  fclose(fids);
  qtmp = transpose(qtmp);
  qtmp  = reshape(qtmp,mx_old,meqn,meth1);
  qsoln = zeros(mx,1,meqn);
  for i=1:mx_old
    for me=1:meqn
      for ii=1:points_per_dir
        v1(1:meth1,1) = phi(ii,:);
        v2(1:meth1,1) = qtmp(i,me,:);
        qsoln((i-1)*points_per_dir+ii,1,me) = transpose(v1)*v2;
      end
    end
  end
  clear qtmp;
  
  qtrack(:,n1+1,:) = qsoln(:,1,:);
  
  %% Aux variables -- aux
  if (maux>0)
    fname(8) = 'a';	      
    fids = fopen(fname,'r');
    tmp = fscanf(fids,'%e',1);
    atmp = fscanf(fids,'%e',[1,inf]);
    fclose(fids);
    atmp = transpose(atmp);
    atmp  = reshape(atmp,mx_old,maux,meth1);
    aux = zeros(mx,1,maux);
    for i=1:mx_old
      for ma=1:maux
        for ii=1:points_per_dir              
          v1(1:meth1,1) = phi(ii,:);
          v2(1:meth1,1) = atmp(i,ma,:);
          aux((i-1)*points_per_dir+ii,1,ma) = transpose(v1)*v2;
        end
      end
    end
  end
  clear atmp;
  
  auxtrack(:,n1+1,:) = aux(:,1,:);
end

[XX,TT]=meshgrid(xc,time);
XX = XX';
TT = TT';
size(XX)
size(TT)

% USER SUPPLIED FUNCTION
plotq1_spacetime;        

disp(' ')
