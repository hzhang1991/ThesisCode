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



figure('Position', [100, 100, 1049, 895]);
clf;
shad = 'flat';
cont = 'off';

x=linspace(0.0,3.0,480);
y=linspace(0.0,1.0,160);
[X,Y]=meshgrid(x,y);
%X=X(:);Y=Y(:);


cx=[node(tnode(:,1),1) node(tnode(:,2),1) node(tnode(:,3),1)];
cx=sum(cx,2)/3.0;
cy=[node(tnode(:,1),2) node(tnode(:,2),2) node(tnode(:,3),2)];
cy=sum(cy,2)/3.0;

[Xq,Yq,Vq] = griddata(cx,cy,qsoln(:,m),X,Y);

x1=Xq(1,:);
y1=Yq(:,1);
ix=find(x1>0.6);
iy=find(y1<=0.2);
[IX,IY]=meshgrid(ix,iy);

Vq(160*(IX(:)-1)+IY(:))=0.0*Vq(IX(:)+IY(:));

colormap(flipud(gray(2048)).^10)
[qx,qy] = gradient(log(Vq),1.0/160.0,1.0/160.0);
qs = sqrt(qx.^2 + qy.^2);
qmesh = qs';
qs(160*(IX(:)-1)+IY(:))=0.0*Vq(IX(:)+IY(:));

h=pcolor(Xq,Yq,qs);
hold on
plot([0.6 3.0],[0.2,0.2],'k','LineWidth',5)
plot([0.6 0.6],[0.0,0.2],'k','LineWidth',5)
set(h, 'EdgeColor', 'none');
%axis on; box off; grid off;
axis('equal');
axis([-0.05 3.05 -0.05 1.05]);
set(gca,'xtick',-4:0.5:4);
set(gca,'ytick',-4:0.2:4);
set(gca,'fontsize',16);
t1 = title(['Density',' at t = ',num2str(time)]); 
set(t1,'fontsize',16);
export_fig -transparent density_schlieren_unstructured.jpg

figure('Position', [100, 100, 1049, 895]);
clf;
%hh=pdeplot(node',[],tnode','xydata',qsoln(:,m),'xystyle',shad,...
%           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','off','mesh','off');
contourf(Xq,Yq,Vq,30);
shading flat;
colormap('jet');
colorbar();
axis on; box on; grid off;
axis('equal');
axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Density',' at t = ',num2str(time)]); 
set(t1,'fontsize',16);
export_fig -transparent density_unstructured.jpg


figure('Position', [100, 100, 1049, 895]);
clf;
H2=pdeplot(node',[],tnode','xydata',transpose(qsoln(:,m)),'xystyle','off',...
           'contour','on','levels',linspace(0.1,4.54,30),'colorbar','off');
for i=1:length(H2)
  set(H2(i),'linewidth',2*0.5);
  set(H2(i),'color','k');
end
axis on; box on; grid off;
axis('equal');
axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['Density',' at t = ',num2str(time)]); 
set(t1,'fontsize',16);
export_fig -transparent density_contour_unstructured.jpg

% 
% figure(3)
% rho=qsoln(:,1);
% u=qsoln(:,2)./qsoln(:,1);
% v=qsoln(:,3)./qsoln(:,1);
% e=qsoln(:,5);
% p=0.4*(e-0.5*rho.*(u.*u+v.*v));
% c=sqrt(1.4*p./rho);
% M=sqrt(u.*u+v.*v)./c;
% colors=linspace(0.5,7.5,60);
% %hh=pdeplot(node',[],tnode','xydata',transpose(M),'xystyle',shad,...
% %           'zdata',M,'zstyle',cont,'colorbar','off','mesh','off');
% H3=pdeplot(node',[],tnode','xydata',transpose(M),'xystyle','off',...
%            'contour','on','levels',colors ,'colorbar','on');
% axis on; box on; grid off;
% axis([-0.1 1.5 -0.5 0.5])
% axis('equal');
% set(gca,'fontsize',16);
% t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
% set(t1,'fontsize',16);
