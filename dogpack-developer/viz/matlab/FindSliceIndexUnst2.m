function [BndyList, Psort, Mid, Left, Right] = FindSliceIndexUnst2( meshdir, outputdir )
%FINDSLICEINDEXUNST2    Find the indices for a given slice through a grid.
%
% This function traverses the grid, and saves the indices for each triangle that
% intercepts the line y=0.
%
% Output is saved in the file, './slice_index.dat'.
%
% Input:
%
%    meshdir - optional directory for where to find the mesh information.
%    Default: './Unstructured_Mesh/mesh_output'.
%
% Output:
%
%    BndyList - indices of triangles that intersect the line y=0.
%
%    Psort    - indices of sorted triangles.
%
% See also: plotq2_unst

% TODO - check if the output information located in outputdir '/slice_index.dat'
% already exists.  If so, just load the data!
%
% fid = fopen([outputdir '/slice_index.dat'],'w');
%

if( nargin < 1 )
  meshdir = ['./Unstructured_Mesh/mesh_output'];
end
if( nargin < 2 )
    outputdir = './Unstructured_Mesh/';
end

fids = fopen([meshdir '/mesh_params.dat'],'r');
if fids==-1
  disp(' ');
  error(['File  ',meshdir,'/mesh_params.dat  not found.']);
end

NumElems      = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumPhysElems  = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumGhostElems = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumNodes      = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumPhysNodes  = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumBndNodes   = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumEdges      = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
fclose(fids);

fids = fopen([meshdir '/mesh_tnode.dat'],'r');
tmp = fscanf(fids,'%d %d %d',[3,inf]);
fclose(fids);
tnode = transpose(tmp);
clear tmp;

fids = fopen([meshdir '/mesh_node.dat'],'r');
tmp = fscanf(fids,'%e',[2,inf]);
fclose(fids);
node = transpose(tmp);
clear tmp;

% This stuff isn't currently being used.
%   fids = fopen([meshdir '/mesh_bnd_node.dat'],'r');
%   tmp = fscanf(fids,'%d',[1,inf]);
%   fclose(fids);
%   bnd_node = transpose(tmp);
%   clear tmp;

%   fids = fopen([meshdir '/mesh_area_prim.dat'],'r');
%   tmp = fscanf(fids,'%e',[1,inf]);
%   fclose(fids);
%   area = transpose(tmp);
%   clear tmp;

%   fids = fopen([meshdir '/mesh_jmat.dat'],'r');
%   Jmat = zeros(NumElems,2,2);
%   for i=1:NumElems
%       tmp = fscanf(fids,'%e',[1,4]);
%       Jmat(i,1,1) = tmp(1);
%       Jmat(i,1,2) = tmp(2);
%       Jmat(i,2,1) = tmp(3);
%       Jmat(i,2,2) = tmp(4);
%   end
%   fclose(fids);
%   clear tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find every triangle that intersects the line y = 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BndyList = [];    
for i=1:NumPhysElems

  % Compute the centroid and location of each vertex (node)
  x = node( tnode(i,:), 1 );   xmid = sum( x ) / 3;
  y = node( tnode(i,:), 2 );   ymid = sum( y ) / 3;

  % compute sgn of the distance from the vertices to the line.  E.g.,
  % determine if the poitns are 'above' or 'below' the line:
  tmp    = sign( y );
  I      = find( abs( y ) < 1e-14 );
  tmp(I) = 0;

  % determine the 'distance' to the line.
  % triangles within '3' of the line contain part of the line on their
  % interior
  d = abs( sum( tmp ) );

  if( d < 3 )
    %    t1=text(xmid,ymid,[num2str(i)]);
    %    set(t1,'color',[0 0 1]);
    %    set(t1,'fontweight','bold');
    BndyList = [BndyList i];
  end

end

% Number of triangles identified as being on the line:
NumTriOnLine = length( BndyList ); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the "midpoint" of each flagged triangle.
%
% While we're at it, compute the 'left' and 'right' values that intersect  the
% boundary of the triangle as well as the line y=0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-10;
Left  = zeros( NumTriOnLine, 2 ); 
Right = zeros( NumTriOnLine, 2 ); 
Mid   = zeros( NumTriOnLine, 2 );
for i=1:NumTriOnLine

  % Compute the centroid and location of each vertex (node)
  x = node( tnode(BndyList(i),:), 1 );   xmid = sum( x ) / 3;
  y = node( tnode(BndyList(i),:), 2 );   ymid = sum( y ) / 3;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Each pair of points on the triangle defines a line.  This line intersects
  % the line y=0 at exactly one point.  Of these three points, only one of them
  % is inside the line.  We will check all three.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pts = [];
  
  % First pair of points:
  x1  = x(1); y1 = y(1);
  x2  = x(2); y2 = y(2);
  tmp = x2 - (x2-x1)/(y2-y1)*y2;
  
  if( min(x1,x2 ) < tmp+tol & tmp-tol < max(x1,x2) )
    pts = [pts; tmp];
  end
  assert( abs(y2-y1) > 1e-13 );
  
  % Second pair of points:
  x1  = x(1); y1 = y(1);
  x2  = x(3); y2 = y(3);
  tmp = x2 - (x2-x1)/(y2-y1)*y2;
  
  if( min(x1,x2 ) < tmp+tol & tmp-tol < max(x1,x2) )
    pts = [pts; tmp];
  end
  assert( abs(y2-y1) > 1e-13 );
  
  % Third pair of points:
  x1  = x(2); y1 = y(2);
  x2  = x(3); y2 = y(3);
  tmp = x2 - (x2-x1)/(y2-y1)*y2;
  
  if( min(x1,x2) < tmp+tol & tmp-tol < max(x1,x2) )
    pts = [pts; tmp];
  end
  assert( abs(y2-y1) > 1e-13 );

  % Sort the two points that we found
  if( ~(length(pts) == 2) )
    % Whoops! you found too many points!
    if( abs( pts(1) - pts(2) ) < tol )
      pts = pts(2:end);
    else
      pts = pts(1:2);
    end
  end
  if( pts(1) > pts(2) )
    pts = pts(2:-1:1);
  end
  
  % Save the points that we found.  These points have already been sorted.
  Left (i,:) = [pts(1); 0.];
  Right(i,:) = [pts(2); 0.];
  Mid  (i,:) = [0.5*( pts(1)+pts(2) ); 0];
  
end
  
% Determine a sorting of 1:NumTriOnLine
% For now, we'll use simple bubble-sort; if need be, this can be sped up.
Psort = 1:NumTriOnLine';
for j=1:NumTriOnLine
  for i=1:(NumTriOnLine-j)
    if( Mid( Psort(i), 1 ) > Mid( Psort(i+1), 1 ) )
      % Swap these two:
      tmp        = Psort(i);
      Psort(i)   = Psort(i+1);
      Psort(i+1) = tmp;
    end
  end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Slice Index
fid = fopen([outputdir '/slice_index.dat'],'w');
for i=1:NumTriOnLine
  fprintf(fid,'%8i %8i\n', BndyList(i), Psort(i) );
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
