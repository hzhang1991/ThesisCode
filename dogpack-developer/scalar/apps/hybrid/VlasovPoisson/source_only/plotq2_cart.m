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

clear qex;
for j=1:my
  for i=1:mx
    rr = sqrt((xc(i,j)-0.25)^2 + (yc(i,j)-0.5)^2);
    if (rr<0.2)
      qex(i,j) = cos(5.0/2.0*pi*rr).^6;
    else
      qex(i,j) = 0.0;
    end
  end
end

figure(1);
clf;
surf(xl,yl,qaug(:,:,1));
colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-0.04 1.04 -0.04 1.04]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 1]);

if (abs(time)<=1.0e-12 || abs(time-1)<=1.0e-12)
  err = reshape(abs(qex-qsoln),mx*my,1);
  err_scale = reshape(abs(qex),mx*my,1);
  
  err_rel = norm(err,2)/norm(err_scale,2);
  
  disp(' ');
  disp([' 2-norm error = ',num2str(err_rel,'%0.5e')]);
  disp(' ');
end
    
figure(2);
clf;
contour(xc,yc,qsoln(:,:,1),11,'k');
axis on; box on; grid off;
axis('equal');
axis([-0.04 1.04 -0.04 1.04]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16); 

figure(3);
clf;
pt=plot(xc(:,round(my/2+1)),qex(:,round(my/2+1)),'r-');
set(pt,'linewidth',1.5);
hold on; 
pz=plot(xc(:,round(my/2+1)),qsoln(:,round(my/2+1),1),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold off; 
axis on; box on; grid off;
axis([0 1 -0.1 1.1]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
yslice = yc(1,round(my/2 + 1));
t1 = title(['q(t,x,',num2str(yslice),') at t = ',...
            num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

figure(1)
