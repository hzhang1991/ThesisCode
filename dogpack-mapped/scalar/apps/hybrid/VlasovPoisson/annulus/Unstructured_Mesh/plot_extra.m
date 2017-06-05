%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EXTRA FUNCTION, SPECIFIED BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

%for i=1:NumNodes
%  xn=node(i,1);
%  yn=node(i,2);
  
%  t1=text(xn,yn,[num2str(i)]);
%  set(t1,'color',[0 0 0]);
%  set(t1,'fontweight','bold');
%  set(t1,'fontsize',16);
%end

for i=1:NumElems
  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;

  t1=text(xmid,ymid,[num2str(i)]);
  set(t1,'color',[0 0 1]);
  set(t1,'fontweight','bold');
  set(t1,'fontsize',12);
end

%for i=1:NumEdges
%  xmid = (edge(i,1)+edge(i,3))/2;
%  ymid = (edge(i,2)+edge(i,4))/2;

%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 1 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',12);
%end

%for i=1:NumPhysElems
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  t1=text(xmid,ymid,[num2str(i)]);
%  set(t1,'color',[0 0 1]);
%  set(t1,'fontweight','bold');
%end

%for j=1:NumGhostElems
%  i = NumPhysElems+j;
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  k = ghost_link(j);
%  
%  t2=text(xmid,ymid,[num2str(k)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%end

%for i=1:NumEdges
%  xmid = (edge(i,1)+edge(i,3))/2;
%  ymid = (edge(i,2)+edge(i,4))/2;
%
%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',14);
%end

%figure(1)
%hold on;
%s=transpose(linspace(0,0.3,20));
%for i=1:NumEdges
%  
%  x1 = edge(i,1);
%  x2 = edge(i,3);
%  y1 = edge(i,2);
%  y2 = edge(i,4);
%  
%  xmid = (x1+x2)/2;
%  ymid = (y1+y2)/2;
%  
%  L = sqrt((x2-x1)^2 + (y2-y1)^2);
%  n1 = (y2-y1)/L;
%  n2 = (x1-x2)/L;
% 
%  xv=xmid+s*n1;
%  yv=ymid+s*n2;
%
%  pvq=plot(xv,yv,'k--');
%  set(pvq,'linewidth',2);
%  
%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',14);
%end
%hold off;
%


%for i=1:NumPhysElems
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  t1=text(xmid,ymid,[num2str(i)]);
%  set(t1,'color',[0 0 1]);
%  set(t1,'fontweight','bold');
%  set(t1,'fontsize',14);
%end

%for i=1:NumPhysNodes
%  t1=text(node(i,1),node(i,2),[num2str(i)]);
%  set(t1,'color',[1 0 1]);
%  set(t1,'fontweight','bold');
%  set(t1,'fontsize',14);
%end

%for i=1:NumEdges
%  xmid = (edge(i,1)+edge(i,3))/2;
%  ymid = (edge(i,2)+edge(i,4))/2;
%
%  t2=text(xmid,ymid,[num2str(i)]);
%  set(t2,'color',[1 0 0]);
%  set(t2,'fontweight','bold');
%  set(t2,'fontsize',14);
%end

figure(1)
