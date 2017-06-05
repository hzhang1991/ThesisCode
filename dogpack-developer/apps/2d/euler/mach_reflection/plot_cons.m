clear
clear functions

format long e

% INFO --------------------------------------------------------------
disp(' ')
disp('GLOBAL CONSERVATION PLOTTER')
disp('    written by James Rossmanith')
disp(' ')
clear;
clear functions;

% meqn and nplot
fids = fopen('output/qhelp.dat','r');
meqn = fscanf(fids,'%d',1);
nplot = fscanf(fids,'%d',1);
fclose(fids);

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





