function plotasynch1(points_per_dir,outputdir_in)

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
%disp(' ');

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

disp(['          melems = ', int2str(mx)]);
disp(' ');

% Grid information
fidg1 = fopen([outputdir,'/mesh_node.dat'],'r');
x = transpose(fscanf(fidg1,'%e',[1 inf]));
mx = length(x)-1;
fclose(fidg1);
dx(1:mx,1) = x(2:(mx+1),1)-x(1:mx);
xc(1:mx,1) = 0.5*(x(2:(mx+1),1)+x(1:mx));
mx_old = mx;
xc_old = xc;
mx = points_per_dir*mx_old;

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
	% solution should be found in file
	%     outputdir/q[n1].dat
    fname = [outputdir,'/',num2str(n1+10000),'.dat'];

    % replace the 1000's digit by the letter q
    fname(length(outputdir)+2) = 'q';

    fids = fopen(fname,'r');
    if fids==-1
      error(['File  ',fname,'  not found.']);
    end

    time = fscanf(fids,'%e',1);
    qtmp = fscanf(fids,'%e',[1,inf]);
    fclose(fids);
    qtmp = transpose(qtmp);
    qcoeffs  = reshape(qtmp,mx_old,meqn,meth1);
    clear qtmp;
    qsoln = zeros(mx,meqn);
    for i=1:mx_old
      for ii=1:points_per_dir
        for me=1:meqn
          v1(1:meth1,1) = phi(ii,:);
          v2(1:meth1,1) = qcoeffs(i,me,:);
          qsoln((i-1)*points_per_dir+ii,me) = transpose(v1)*v2;
        end
        xc((i-1)*points_per_dir+ii,1) = x(i) + (ii-0.5)*dx(i)/(points_per_dir);
      end
    end
    
    %% Aux variables -- aux
    if (maux>0)
      fname(length(outputdir)+2) = 'a';
      fids = fopen(fname,'r');
      if fids==-1
        error(['File  ',fname,'  not found.']);
      end
      time = fscanf(fids,'%e',1);
      atmp = fscanf(fids,'%e',[1,inf]);
      fclose(fids);
      atmp = transpose(atmp);
      acoeffs  = reshape(atmp,mx_old,maux,meth1);
      clear atmp;
      aux = zeros(mx,maux);
      for i=1:mx_old
        for ma=1:maux
          for ii=1:points_per_dir              
            v1(1:meth1,1) = phi(ii,:);
            v2(1:meth1,1) = acoeffs(i,ma,:);
            aux((i-1)*points_per_dir+ii,ma) = transpose(v1)*v2;
          end
        end
      end
    end
      
    % USER SUPPLIED FUNCTION
    plotq1;        
  end

end
disp(' ')
