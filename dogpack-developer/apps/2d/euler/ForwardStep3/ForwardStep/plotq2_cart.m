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

% gas constant
fids  = fopen([outputdir,'/eulerhelp.dat'],'r');
if fids==-1
  disp(['File  ',outputdir,'/eulerhelp.dat  not found.']);
  disp('Setting gamma = 1.4');
  gamma_gas = 1.4;
else
  gamma_gas  = fscanf(fids,'%e',1);
  fclose(fids);
end

figure(1);
clf;
yrbcolormap
rho=qsoln(:,:,m)/1.4;
pcolor(xc,yc,rho);
%contourf(xc, yc, rho, linspace(0.1,4.54,30), '-k' );
shading flat;
axis on; box on; grid off;
axis('equal');
axis([-0.05 3.05 -0.05 1.05]);
set(gca,'xtick',-4:1.0:14);
set(gca,'ytick',-4:1.0:14);
set(gca,'fontsize',16);
time
t1 = title(['Density'])
%t1 = title(['Density at t = ',num2str(time)']); 
set(t1,'fontsize',16);
colorbar;
%caxis([0.8,10]);
caxis auto;
export_fig(sprintf('pdfs2/density_color%02d.png',n1),'-transparent')% disp(['density_color',num2str(nf),'.pdf'])


% % Grey scale figure plot that you see in a lot of papers.
figure(2);
clf;
contourf(xc, yc, rho, 40, '-k' );
colorbar()
axis on; box on; grid off;
axis('equal');
axis([-0.05 3.05 -0.05 1.05]);
set(gca,'xtick',-4:1.0:14);
set(gca,'ytick',-4:1.0:14);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [DOGPACK]']); 
set(t1,'fontsize',16);

figure(3)
clf;
rho=qsoln(:,:,1);
e=qsoln(:,:,5);
u1=qsoln(:,:,2)./qsoln(:,:,1);
u2=qsoln(:,:,3)./qsoln(:,:,1);
u3=qsoln(:,:,4)./qsoln(:,:,1);
P=0.4*(e-0.5*rho.*(u1.*u1+u2.*u2+u3.*u3));
P(xc<=1 & yc<=6)=0.0;
c=sqrt(1.4*P./rho);
M=sqrt(u1.*u1+u2.*u2)./c;
%contour(xc, yc, P, linspace(0.091, 38, 50), '-k' );
contour(xc, yc, rho, linspace(0.1,4.54,30), '-k' );
colorbar()
axis on; box on; grid off;
axis('equal');
axis([-0.05 3.05 -0.05 1.05]);
set(gca,'xtick',-4:1.0:14);
set(gca,'ytick',-4:1.0:14);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time)]); 
set(t1,'fontsize',16);
export_fig -transparent density_contour.pdf
% n1 = frame number
% fname = strcat( strcat( 'density', num2str(n1, '%02d' ) ), '.jpg' );
% print(1, '-djpeg', fname );
% fname = strcat( strcat( 'density-contour', num2str(n1, '%02d' ) ), '.eps' );
% print(3, '-deps', fname  );
