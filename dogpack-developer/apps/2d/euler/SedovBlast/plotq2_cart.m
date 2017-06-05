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

% ------------------------------------------------------------- %
% Density (across the entire domain)
% ------------------------------------------------------------- %
figure(1);
clf;
v=linspace(0,6,40);
contourf(xl,yl,qaug(:,:,m),v);
shading flat;
yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-1.1  1.1  -1.1  1.1 ]);
set(gca,'xtick',-4:0.5:4);
set(gca,'ytick',-4:0.5:4);
set(gca,'fontsize',16);
t1 = title(['\rho(t = ', num2str(time,'%2.2e'), ', x)     [DoGPack]']); 
set(t1,'fontsize',16);
colorbar;
caxis([0,6]);

% Construct scatter plot of current solution - save as qscat
rscat = zeros( mx*my, 1    );
qscat = zeros( mx*my, meqn );
k = 1;
for i=1:mx
for j=1:my
    qscat(k, : ) = qsoln(i,j,:);
    rscat(k    ) = sqrt( xc(i,j)^2 + yc(i,j)^2 );
    k = k+1;
end
end

% ------------------------------------------------------------- %
% Density (scatter plot for radially symmetric density)
% ------------------------------------------------------------- %
figure(2)
clf
pz = plot( rscat(:), qscat(:,1), 'bo' );
set(pz, 'linewidth',  1);
set(pz, 'markersize', 6);
t1 = title(['Density (scatter plot)']);
set(t1, 'fontsize', 16);
axis on; box on; grid off;

% Add in the "exact" solution.  This data was pulled from Sedov's original
% text.
r   = [0 .1 .15 .2 .25 .3 .3654 .4222 .4748 .518 .5754 .639 .6894 .7274 .7629 .8094 .8442 .8725 .9096 .9295 .9476 .9644 .9802 .998 1];
den = [ 0 0 0 0.0001 .0003 .0008 .0021 .0044 .0079 .0123 .0208 .0362 .0545 .0718 .0975 .1414 .1892 .2427 .3451 .4234 .5164 .6285 .7653 .9973 1];

% Additional parameters
rho = 1.0;
d0  = (1.4+1.)/(1.4-1.)*rho;
r0  = 1.;
dd = den*d0;
rr = r*r0;
rr = [rr r0 1.1];
dd = [dd 1 1 ];

hold on;
pz2 = plot( rr, dd, '-r' );
set( pz2, 'linewidth', 1 );
hold off;
axis([0 1.1 0 6.5]);

figure(1)

% n1 = frame number
%fname = strcat( strcat( 'density', num2str(n1, '%02d' ) ), '.jpg' );
%print(1, '-djpeg', fname );
%fname = strcat( strcat( 'density-contour', num2str(n1, '%02d' ) ), '.eps' );
%print(3, '-deps', fname  );

