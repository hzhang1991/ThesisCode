%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%              mx:  number of points
%%           nplot:  number of frames
%%    [xlow,xhigh]:  min/max values of grid
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy           
%%
%%   Grid information:
%%         [XX,TT]: spacetime mesh (both XX and TT are size = (mx,nplot+1))
%%
%%   Solution information:
%%          qtrack:  solution sampled on spacetime mesh, size = (mx,nplot+1,meqn)
%%        auxtrack:  aux components sampled on mesh, size = (mx,nplot+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmin = min(min(XX));
xmax = max(max(XX));
tmin = min(min(TT));
tmax = max(max(TT));

figure(1);
surf(XX,TT,qtrack(:,:,m));
shading flat;
view([-35 80]);
set(gca,'xlim',[xmin xmax]);
set(gca,'ylim',[tmin tmax]);
colormap('jet');
set(gca,'fontsize',16);
t1=title(['Spacetime plot of q(',num2str(m),')']);
set(t1,'fontsize',16);
