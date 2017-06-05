function plotdog2_cons(outputdir_in)

global outputdir
if(nargin<1)
  outputdir='output';
else
  outputdir=outputdir_in;
end

format long e

% INFO --------------------------------------------------------------
disp(' ')
disp('2D DoGPack CONSERVATION PLOTTER')
disp('    written by the DoGPack team')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse the file: QHELP.DAT
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
else
  error(['Incorrect GridType, must be Cartesian or Unstructured. GridType = ',GridType]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(['        GridType = ',GridType]);
disp(['       outputdir = ',outputdir]);
disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CARTESIAN conservation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unstructured conservation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (GridType=='Unstructured')  

  meqn          = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  maux          = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  nplot         = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  meth1         = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  datafmt       = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  %
  NumElems      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumPhysElems  = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumGhostElems = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumNodes      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumPhysNodes  = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumBndNodes   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumEdges      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);  
  fclose(fids);
 
else
  error([' Unsupported value: GridType = ',GridType]);
end

m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
	' ) ? ']);
disp(' ')
if isempty(m)
  m=1;
end

% INFO --------------------------------------------------------------

fid = fopen([outputdir '/conservation.dat']);
consv = fscanf(fid,'%e',[meqn+1 inf]);
status = fclose(fid);

t = consv(1,:)';
qc = consv(m+1,:)';
nt = length(t);

%dqc = qc - qc(1);
for i=1:nt
  dqc(i) = qc(i)-qc(1);
end

figure(2)
clf
plot(t,dqc,'b-');
set(gca,'fontsize',16);
t1=title([ 'Conservation of q(',num2str(m),') vs. time']);
set(t1,'fontsize',16);
axis on; box on; grid off;
hold off
