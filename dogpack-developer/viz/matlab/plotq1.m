%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%              mx:  number of points
%%    [xlow,xhigh]:  min/max values of grid
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%              xc: grid points (cell centers), size = (mx,my)
%%
%%   Solution information:
%%           qsoln:  solution sampled on mesh, size = (mx,meqn)
%%             aux:  aux components sampled on mesh, size = (mx,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qlow  = min(qsoln(:,m));
qhigh = max(qsoln(:,m));
qeps  = 0.015*(qhigh-qlow);

figure(1);
clf;
plot(xc,qsoln(:,m));
axis on; box on; grid off;
axis([xlow xhigh qlow-qeps qhigh+qeps]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(1),') at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
