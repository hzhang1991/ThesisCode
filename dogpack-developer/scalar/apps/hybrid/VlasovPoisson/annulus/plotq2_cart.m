%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                  mx, my:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh]:  min/max values of grid
%%                    meqn:  number of equations
%%                    maux:  number of aux components
%%                   meth1:  spatial order of accuracy
%%
%%   Grid information:
%%       (xc,yc): grid points (cell centers), size = (mx,my)
%%       (xl,yl): grid points (lower left cell corners), size = (mx+1,my+1)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,maux)
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,meqn)
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% background maxwellian:
%
% units with meaning:
N = 1.0e13; % 1/m
L = 0.1; % m
T = 5.605424055e-9; % s
E0 = 1.809512620e4; % N/C
Phi0 = 1.809512620e3; % Nm/C

% was this already here?
temp = 5.526350206e-4;


%max( max( qsoln ) )
%min( min( qsoln ) )

figure(1);
clf;
pcolor(xl,yl,qaug(:,:,1));
colormap('jet');
axis on; %box on; grid off;
grid off;
axis('equal');
axis([-0.21 0.21 -0.21 0.21]);
%set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['f(t,vx,vy) at t = ',num2str(T*time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 300]);
colorbar;
shading flat;

   
% plot the difference between this solution and a maxwellian
rho0 = 2.830722109641852e+02 * (2*pi*temp);
maxwellian = rho0 * exp( -(xc.^2+yc.^2)/(2*temp) ) / (2*pi*temp);

figure(2);
clf;
surf(xc,yc, maxwellian - qsoln );
colormap('jet');
axis on; %box on; grid off;
grid off;
axis('equal');
axis([-0.21 0.21 -0.21 0.21]);
%set(gca,'xtick',-2:0.25:2);
%set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['f(t,vx,vy) at t = ',num2str(T*time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([-80 80 ] );
caxis([-28 28 ] );
%caxis auto;
colorbar;
shading flat;

%   figure(2);
%   clf;
%   contour(xc,yc,qsoln(:,:,1),11,'k');
%   axis on; box on; grid off;
%   axis('equal');
%   axis([-0.04 1.04 -0.04 1.04]);
%   set(gca,'xtick',-2:0.25:2);
%   set(gca,'ytick',-2:0.25:2);
%   set(gca,'fontsize',16);
%   t1 = title(['f(t,x,y) at t = ',num2str(T*time),'     [DoGPack]']); 
%   set(t1,'fontsize',16); 

figure(3);
clf;
%pt=plot(xc(:,round(my/2+1)),qex(:,round(my/2+1)),'r-');
%set(pt,'linewidth',1.5);
%hold on; 
pz=plot(xc(:,round(my/2+1)),qsoln(:,round(my/2+1),1),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold off; 
axis on; box on; grid off;
axis([-0.21 0.21 -0.1 290]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-0.2:0.05:0.2);
set(gca,'ytick',0:30:300);
set(gca,'fontsize',16);
yslice = yc(1,round(my/2 + 1));
t1 = title(['f(t,vx,',num2str(yslice),') at t = ',...
            num2str(T*time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(1)

  % print the pretty pictures!
  % frame number = n1
% folder_name = strcat( outputdir, '/photos/' );
% fname = strcat( strcat('sheath-2d-phase-diff', num2str(n1, '%03d') ), '.eps' );
% print(2,'-depsc', strcat(folder_name, fname ) );
