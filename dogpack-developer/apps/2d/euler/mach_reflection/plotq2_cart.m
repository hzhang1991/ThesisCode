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
  error(['File  ',outputdir,'/eulerhelp.dat  not found.']);
end
gamma_gas  = fscanf(fids,'%e',1);
fclose(fids);

rho=qsoln(:,:,m);
levels=linspace(1.5, 22.7,30);

figure('Position', [100, 100, 1049, 349.67]);
clf;
%linspace(1.1, 22,30)
contour(xc,yc,rho,levels,'-k');
   set(gca,'xtick',-4:0.5:4);
   set(gca,'ytick',-4:0.5:4);
 t1 = title(['Density at t = ',num2str(time)]); 
set(t1,'fontsize',16);
%export_fig -transparent contour_p3.jpg
export_fig -transparent contour_p3.pdf

figure('Position', [100, 100, 349.67, 349.67]);
clf;
%linspace(1.1, 22,30)
contour(xc,yc,rho,levels,'-k');
axis('equal');
axis([2.25 3.0 -0.05 0.75]);
   set(gca,'xtick',-4:0.5:4);
   set(gca,'ytick',-4:0.5:4);
t1 = title(['Density at t = ',num2str(time)]); 
set(t1,'fontsize',16);

%export_fig -transparent contour_p3_zoom.jpg
export_fig -transparent contour_p3_zoom.pdf

