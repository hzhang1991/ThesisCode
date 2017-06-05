%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%                NumElems:  number of elements
%%            NumPhysElems:  number of elements (excluding ghost elements)
%%           NumGhostElems:  number of ghost elements
%%                NumNodes:  number of nodes
%%            NumPhysNodes:  number of nodes (excluding exterior to domain)
%%             NumBndNodes:  number of nodes on boundary
%%                NumEdges:  number of edges
%% [xlow,xhigh,ylow,yhigh]:  bounding box
%%                    node:  list of nodes in mesh
%%                   tnode:  nodes attached to which element
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (NumPhysElem,meqn)
%%           aux:  aux components sampled on mesh, size = (NumPhysElem,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%xstar = -0.25*cos(2.0*pi*time) + 0.50;
%ystar =  0.25*sin(2.0*pi*time) + 0.50;
xstar = 0.0;
ystar = 0.0;

for i=1:(NumPhysElems*points_per_dir^2)
  xm = xmid(i);
  ym = ymid(i);

  rscat(i,1) = sqrt((xm-xstar)^2 + (ym-ystar)^2);
% if (rscat(i,1)<0.2)
%   qex(i,1) = cos(5.0/2.0*pi*rscat(i,1)).^6;
% else
%   qex(i,1) = 0.0;
% end
    % qex(i,1) = (0.75+0.25*cos(2*pi*time)) * (0.5+0.5*cos(pi*rscat(i,1))); 
    qex(i,1) = (0.5+0.5*cos(pi*rscat(i,1))); 
end

[rscat,Itmp]=sort(rscat);
qex(:,1) = qex(Itmp,1);


figure(1);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,m),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','off','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-1.01 1.01 -1.01 1.01]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 1]);
hold off;

err = norm(qex-qsoln(Itmp,1),2)/norm(qex,2);
disp(' ');
disp(['  Error = ',num2str(err,'%0.5e')]);
disp(' ');

figure(2)
clf;
H2=pdeplot(node',[],tnode','xydata',qsoln(:,m),'xystyle','off',...
           'contour','on','levels',12,'colorbar','off');
for i=1:length(H2)
  set(H2(i),'linewidth',2*0.5);
  set(H2(i),'color','k');
end
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(3);
clf;
pz=plot(rscat,qsoln(Itmp,m),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold on;
pt=plot(rscat,qex(:,m),'r-');
set(pt,'linewidth',2);
hold off; 
axis on; box on; grid off;
axis([0 1.0 -0.1 1.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(t,r) at t = ',...
            num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(1)

