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

figure(1);
clf;
shad = 'flat';
cont = 'off';
hh=pdeplot(node',[],tnode','xydata',qsoln(:,m),'xystyle',shad,...
           'zdata',qsoln(:,m),'zstyle',cont,'colorbar','off','mesh','off');
colormap('jet');
colorbar();
axis on; box on; grid off;
axis('equal');
axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

% figure(2)
% clf;
% H2=pdeplot(node',[],tnode','xydata',transpose(qsoln(:,m)),'xystyle','off',...
%            'contour','on','levels',12,'colorbar','off');
% for i=1:length(H2)
%   set(H2(i),'linewidth',2*0.5);
%   set(H2(i),'color','k');
% end
% axis on; box on; grid off;
% axis('equal');
% axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
% set(gca,'fontsize',16);
% t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
% set(t1,'fontsize',16);
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
