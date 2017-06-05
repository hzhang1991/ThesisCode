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

% TODO - use these units!
% units with meaning:
N = 1.0e13; % 1/m
L = 0.1; % m
T = 5.605424055e-9; % s
E0 = 1.809512620e4; % N/C
Phi0 = 1.809512620e3; % Nm/C

disp([['Largest value of q seen = ', num2str(max(max(qsoln)))]]);

xstar = 0.;
ystar = 0.;

for i=1:(NumPhysElems*points_per_dir^2)
  xm = xmid(i);
  ym = ymid(i);
  rscat(i,1) = sqrt((xm-xstar)^2 + (ym-ystar)^2);
end

[rscat,Itmp]=sort(rscat);

figure(1);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(L*node',[],tnode','xydata',qsoln(:,m),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','off','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis('equal');
%axis([-0.01 1.01 -0.01 1.01]);
%set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['rho(t,x,y) at t = ',num2str(T*time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 1]);
hold off;

%err = norm(qex-qsoln(Itmp,1),2)/norm(qex,2);
%disp(' ');
%disp(['  Error = ',num2str(err,'%0.5e')]);
%disp(' ');

figure(2)
clf;
H2=pdeplot(L*node',[],tnode','xydata',qsoln(:,m),'xystyle','off',...
           'contour','on','levels',12,'colorbar','off');
for i=1:length(H2)
  set(H2(i),'linewidth',2*0.5);
  set(H2(i),'color','k');
end
axis on; box on; grid off;
axis('equal');
%axis([-0.01 1.01 -0.01 1.01]);
%set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y) at t = ',num2str(T*time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(3);
clf;
pz=plot(L*rscat, N*qsoln(Itmp,m),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold on;
pt=plot(L*rscat, N*ones(size(rscat)), '--k');
set(pt,'linewidth',2);
hold off; 
axis on; box on; grid off;
%axis([1 10.1 -0.1 1.25e13]);
%set(gca,'plotboxaspectratio',[1.5 1 1]);
%set(gca,'xtick',0:0.02:0.1);
%set(gca,'ytick',0:0.25e13:1.25e13);
set(gca,'fontsize',16);
t1 = title(['n_e(t,r) at t = ',...
            num2str(T*time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(1)

  % print the pretty pictures!
  % frame number = n1
% folder_name = strcat( outputdir, '/photos/' );
% fname = strcat( strcat('sheath-2d-phase', num2str(n1, '%03d') ), '.eps' );
% print(1,'-depsc', strcat(folder_name, fname ) );

% fname = strcat( strcat('sheath-2d-radial', num2str(n1, '%03d') ), '.eps' );
% print(3,'-depsc', strcat( folder_name, fname ) );

