%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%              mx:  number of points
%%    [xlow,xhigh]:  min/max values of grid
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%              xc: grid points (cell centers), size = (mx,my)
%%
%%   Solution information:
%%           qsoln:  solution sampled on mesh, size = (mx,meqn)
%%             aux:  aux components sampled on mesh, size = (mx,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gas constant
%% ----
%% HACK
%% ----
gamma = 1.4;
%% ----


 fname = ['output_ref','/',num2str(20+10000),'.dat'];
      % replace the 1000's digit by the letter q
 fname(length('output_ref')+2) = 'q';
fids  = fopen(['output_ref','/qhelp.dat'],'r');
 
 if fids==-1
          error(['File  ',fname,'  not found.']);
      end
      
 if( fids==-1 )
      error(['File  ','output_ref','/qhelp.dat  not found.']);
 end
 ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
 if( ndims~=1 )
      error(['Incorrect dimension, ndims must be 1. ndims = ',num2str(ndims)]);
 end
 meqn1    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
 maux1    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
 nplot1   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  %
 mx1      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
 xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
 xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
       %% Conserved variables -- qsoln
      fids = fopen(fname,'r');
      if fids==-1
          error(['File  ',fname,'  not found.']);
      end
 

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Read in the data from the output folder
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %% Conserved variables -- qsoln
      time  = fscanf(fids,'%e', 1);
      qtmp  = fscanf(fids,'%e', [1,inf]);
      %qsoln = reshape( qtmp, mx1, meqn1 );

   dx1 = (xhigh-xlow)/mx1;                                 % cell spacing

   xc1 = transpose(linspace(xlow+dx1/2, xhigh-dx1/2, mx1));  % cell centers 
   mx1
   meqn1
   size(qtmp)
      qsolnref = reshape( qtmp, mx1, meqn1 );
      clear qtmp;     

figure(1);
clf;
pz=plot(xc,qsoln(:,1),'b.');
hold on
plot(xc1,qsolnref(:,1),'r-');
legend('DG','Reference','Location','SouthWest')
set(pz,'linewidth',1);
set(pz,'markersize',8);
hold off;
axis on; box on; grid on;
%axis([-5 5 0 5]);
%axis([0.45 2.5 2.75 5.0])
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-5:2.5:5);
set(gca,'ytick',0:0.5:5);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time)]); 
set(t1,'fontsize',16);
export_fig -transparent big_alpha.pdf
%export_fig -transparent small_alpha_averages.pdf
%export_fig -transparent small_alpha.pdf
% figure(2);
% clf;
% press = (gamma-1).*(qsoln(:,5)-0.5*(qsoln(:,2).^2 + ...
%                                     qsoln(:,3).^2 + qsoln(:, ...
%                                                   4).^2)./qsoln(:,1));
% pz=plot(xc,press,'bo');
% set(pz,'linewidth',1);
% set(pz,'markersize',8);
% hold off;
% axis on; box on; grid off;
% %axis([-5 5 0 12]);
% set(gca,'plotboxaspectratio',[1.5 1 1]);
% set(gca,'xtick',-5:2.5:5);
% set(gca,'ytick',0:2:12);
% set(gca,'fontsize',16);
% t1 = title(['Pressure at t = ',num2str(time),'     [DoGPack]']); 
% set(t1,'fontsize',16);
% 
% figure(3);
% clf;
% pz=plot(xc,qsoln(:,2)./qsoln(:,1),'bo');
% set(pz,'markersize',8);
% set(pz,'linewidth',1);
% hold off;
% axis on; box on; grid off;
% %axis([-5 5 -0.5 3]);
% set(gca,'plotboxaspectratio',[1.5 1 1]);
% set(gca,'xtick',-5:2.5:5);
% set(gca,'ytick',0:1:3);
% set(gca,'fontsize',16);
% t1 = title(['u^1(x,t) at t = ',num2str(time),'     [DoGPack]']); 
% set(t1,'fontsize',16);
%     
% figure(4);
% clf;
% pz=plot(xc,qsoln(:,1),'bo');
% set(pz,'linewidth',1);
% set(pz,'markersize',8);
% hold off;
% axis on; box on; grid off;
% %axis([1.0 3.2 2.5 5])
% set(gca,'plotboxaspectratio',[1.5 1 1]);
% set(gca,'xtick',-5:2.5:5);
% set(gca,'ytick',0:0.5:5)
% set(gca,'fontsize',16);
% t1 = title(['Density at t = ',num2str(time),'     [DoGPack]']); 
% set(t1,'fontsize',16);
