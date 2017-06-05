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

% Initial conditions (for Landau damping) are
%   f0 = 0.5/pi * exp( -0.5*(vx*vx+vy*vy) ) * ( 1.0 + alpha*cos(kx*x)*cos(ky*y) );
%
% Exact density (at initial conditions) are given by
%
%     rho0 = 1.0 + alpha*cos( kx*x )*cos( ky*y );
%    rhou0 = 0.0;

% Parameter used for setting up this problem
alpha = 0.01;

% Compare the integrated density of f to the exact initial conditions for rho.
if( n1 == 0 )

    mtmp = NumPhysElems*points_per_dir^2;
    rho0 = zeros(mtmp,1);
    for i=1:mtmp
      xm = xmid(i);
      ym = ymid(i);
      rho0(i,1) = 1.0+alpha*cos( 0.5*xm )*cos( 0.5*ym );
    end

    err = norm(rho0-qsoln(:,1),2)/norm(rho0,2);
    disp(' ');
    disp(['  Error in initial conditions of the density = ',num2str(err,'%0.5e')]);
    disp(' ');
end

%max( qsoln(:,1) )
%min( qsoln(:,1) )

% Plot of whatever the user requested (component passed in to plotdog4)
figure(1);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,m),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','on','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis auto;
axis('equal');
axis([-0.01 12.6 -6.31 6.31]);
set(gca,'xtick',-6:1.00:13);
set(gca,'ytick',-6:1.00:13);
set(gca,'fontsize',16);
%t1 = title(['rho(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
t1 = title(['q', num2str(m), ' at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
%caxis([0 1]);
hold off;


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
axis([-0.01 12.6 -6.31 6.31]);
set(gca,'xtick',-6:1.00:13);
set(gca,'ytick',-6:1.00:13);
set(gca,'fontsize',16);
%t1 = title(['rho(t,x,y) at t = ', num2str(time), '     [DoGPack]']); 
t1 = title(['q', num2str(m), ' at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

% Plot of the log of whateve component was requested
figure(4);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata', log( abs(qsoln(:,m) ) ),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','off','mesh','off');
colormap('jet');
axis on; box on; grid off;
axis auto;
axis('equal');
axis([-0.01 12.6 -6.31 6.31]);
set(gca,'xtick',-6:1.00:13);
set(gca,'ytick',-6:1.00:13);
set(gca,'fontsize',16);
%t1 = title(['q(t,x,y) at t = ', num2str(time), '     [DoGPack]']); 
t1 = title(['log( q', num2str(m), ' ) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
hold off;

%   figure(5);
%   clf;
%   shad = 'flat';
%   cont = 'off';
%   hh=pdeplot(node',[],tnode','xydata',rho0,'xystyle',shad,...
%              'zdata',rho0,'zstyle',cont,'colorbar','on','mesh','off');
%   colormap('jet');
%   axis on; box on; grid off;
%   axis auto;
%   axis('equal');
%   axis([-0.01 12.6 -6.31 6.31]);
%   set(gca,'xtick',-6:1.00:13);
%   set(gca,'ytick',-6:1.00:13);
%   set(gca,'fontsize',16);
%   t1 = title(['rho(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
%   set(t1,'fontsize',16);
%   %caxis([0 1]);
%   hold off;

% print the pretty pictures!
%   Recall : frame number = n1
folder_name = strcat( outputdir, '/photos/' );
if( ~exist(folder_name, 'dir' ) )
    mkdir( folder_name );
end
fname = strcat( strcat('beam-density', num2str(n1, '%03d') ), '.eps' );
print(1,'-depsc', strcat(folder_name, fname ) );

