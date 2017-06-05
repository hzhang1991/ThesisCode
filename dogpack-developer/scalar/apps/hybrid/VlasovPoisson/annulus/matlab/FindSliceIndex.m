function BndyList = FindSliceIndex( )
%FINDSLICEINDEX
%
% This function traverses the grid, and saves the indices for each triangle that
% intercepts the line y=0.
%
% Output is saved in the file, './slice_index.dat'.
%
% Other curves can also be accomodated.

outputdir = '.';
meshdir = ['../Unstructured_Mesh/mesh_output'];

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

fids = fopen([meshdir '/mesh_bnd_node.dat'],'r');
tmp = fscanf(fids,'%d',[1,inf]);
fclose(fids);
bnd_node = transpose(tmp);
clear tmp;

fids = fopen([meshdir '/mesh_area_prim.dat'],'r');
tmp = fscanf(fids,'%e',[1,inf]);
fclose(fids);
area = transpose(tmp);
clear tmp;

fids = fopen([meshdir '/mesh_jmat.dat'],'r');
Jmat = zeros(NumElems,2,2);
for i=1:NumElems
    tmp = fscanf(fids,'%e',[1,4]);
    Jmat(i,1,1) = tmp(1);
    Jmat(i,1,2) = tmp(2);
    Jmat(i,2,1) = tmp(3);
    Jmat(i,2,2) = tmp(4);
end
fclose(fids);
clear tmp;


% Find every triangle that intersects the line y = 0.
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
Ntri = length( BndyList );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Slice Index
fid = fopen([outputdir '/slice_index.dat'],'w');
for i=1:Ntri
    fprintf(fid,'%8i\n', BndyList(i));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
