function plotdog1_cons(outputdir_in)
%PLOTDOG1_CONS     1D Conservation plotter.
%
% Read in and plot the contents of 'outputdir/conservation.dat'
%
% See also: plotdog1, plotdog2_cons

global outputdir
outputdir='output';

if(nargin>1)
  outputdir=outputdir_in;
elseif(isempty(outputdir))
  outputdir='output';
end

format long e

% INFO --------------------------------------------------------------
disp(' ')
disp('1D DoGPack CONSERVATION PLOTTER')
disp('    written by the DoGPack team.')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse the file: QHELP.DAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fids  = fopen([outputdir,'/qhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/qhelp.dat  not found.']);
end
ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
if (ndims~=1)
  error(['Incorrect dimension, ndims must be 1. ndims = ',num2str(ndims)]);
end
GridType = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
GridTmp = GridType(1:9);
if (GridTmp=='Cartesian')
  GridType='Cartesian   ';
else
  error(['Incorrect GridType, in 1D it must be Cartesian. GridType = ',GridType]);
end

meqn    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
maux    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
nplot   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
meth1   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
datafmt = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
%
mx      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
fclose(fids);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(['        GridType = ',GridType]);
disp(['       outputdir = ',outputdir]);
disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
	' ) ? ']);
disp(' ')
if isempty(m)
  m=1;
end

% INFO --------------------------------------------------------------

fid = fopen('output/conservation.dat');
consv = fscanf(fid,'%e',[meqn+1 inf]);
status = fclose(fid);

t = consv(1,:)';
qc = consv(m+1,:)';
nt = length(t);

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
