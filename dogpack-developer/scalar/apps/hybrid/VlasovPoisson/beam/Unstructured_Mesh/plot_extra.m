%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EXTRA FUNCTION, SPECIFIED BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

  % Add in the index for each node (vertex)
% for i=1:NumNodes
%   xn=node(i,1);
%   yn=node(i,2);
%     
%   t1=text(xn,yn,[num2str(i)]);
%   set(t1,'color',[0 0 1]);
%   set(t1,'fontweight','bold');
%   set(t1,'fontsize',14);
% end

 % Plot the centroid of each element
%for i=1:NumElems
%  xmid = (node(tnode(i,1),1)+node(tnode(i,2),1)+node(tnode(i,3),1))/3;
%  ymid = (node(tnode(i,1),2)+node(tnode(i,2),2)+node(tnode(i,3),2))/3;
%
%  t1=text(xmid,ymid,[num2str(i)]);
%  set(t1,'color',[0 0 1]);
%  set(t1,'fontweight','bold');
%end

addpath('../matlab/');
BndyList = FindSliceIndex();
rmpath('../matlab/');

% plot the intersected triangles:
% Can perform a contour fit on curve_func here if we really want ...
hold on;
xline = linspace( xmin, xmax );
p1 = plot( xline, 0.0*xline, '--k');
set(p1, 'LineWidth', 2 );
hold off;
pause;

%% -------------------------------------------------- %%
% search for a single point on the curve:
%% -------------------------------------------------- %%


%   xpt = 0.0;
%   pt  = [xpt, -0.5];

%   % plot the point we're searching for
%   hold on;
%   plot( [pt(1) pt(1)], [pt(2) pt(2)], 'ko', 'LineWidth', 3 );
%   hold off;

%   fprintf(1,'Searching for a single point (%f,%f) in the grid\n', pt(1),pt(2) );
%   NumBndTriangles = length( BndyList );
%   for i=1:NumBndTriangles

%       % Compute the centroid and location of each vertex (node)
%       x = node( tnode(BndyList(i),:), 1 );   xmid = sum( x ) / 3;
%       y = node( tnode(BndyList(i),:), 2 );   ymid = sum( y ) / 3;

%       dx = pt(1) - xmid;  dy = pt(2) - ymid;
%       A  = [ [x(2)-x(1), x(3)-x(1)]; [y(2)-y(1), y(3)-y(1)] ];
%       spts = A \ [dx,dy]';

%       xi = spts(1);  eta = spts(2);

%       % check if this point is inside the triangle:
%       one_third = 1./3.;
%       if( xi > -one_third & eta > -one_third & eta+xi < one_third )
%           fprintf(1,'    found the pt (%f,%f) in triangle number %d\n', ...
%           pt(1), pt(2), BndyList(i));
%       end
%   end

figure(1)
