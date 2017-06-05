%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                         mx, my, mz:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh,zlow,zhigh]:  min/max values of grid
%%                               meqn:  number of equations
%%                               maux:  number of aux components
%%                              meth1:  spatial order of accuracy
%%
%%   Grid information:
%%       (xc,yc,zc): grid points (cell centers), size = (mx,my,mz)
%%       (xl,yl,zl): grid points (lower left cell corners), size = (mx+1,my+1,mz+1)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,mz,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,mz,maux)
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,mz+1,meqn)
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,mz+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
clf;

figure(1);