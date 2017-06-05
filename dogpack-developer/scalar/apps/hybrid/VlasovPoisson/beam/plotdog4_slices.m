% Script written to plot slices of the solution.  This script will read in the
% indices from matlab/slice_index.dat, and plot stuff along the line y=0.

clear

% Parameters that would normally be passed in if this were a function:
global outputdir
outputdir='output';
points_per_dir = 8;   % MUST EQUAL 1
point_type     = 1;  

format long e;

%% ------------------------------------------------------------------------- %%
% This application has a bastardized version of 'unstructured' vs.
% 'structured'.  In order to make this routine work, The paramter that is set
% in qhelp is 'unstructured', because it needs access to a mesh in order to
% perform its plotting duties.
%% ------------------------------------------------------------------------- %%
fids  = fopen([outputdir,'/qhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/qhelp.dat  not found.']);
end
ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
if (ndims~=4)
    error(['Incorrect dimension, ndims must be 4. ndims = ',num2str(ndims)]);
end
GridType = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
GridType = 'Unstructured';
%% ------------------------------------------------------------------------- %%

%% ------------------------------------------------------------------------- %%
% Friendly welcome message:
%% ------------------------------------------------------------------------- %%
disp(' ');
disp(['        GridType = ',GridType]);
disp(['  points_per_dir = ',num2str(points_per_dir)]);
disp(['      point_type = ',num2str(point_type)]);
disp(['       outputdir = ',outputdir]);
disp(' ');
%% ------------------------------------------------------------------------- %%

if (GridType=='Unstructured')  

  % Read in the "Unstructured" data:
  meqn          = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  maux          = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  nplot         = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  meth1         = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  datafmt       = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  %
  NumElems      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumPhysElems  = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumGhostElems = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumNodes      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumPhysNodes  = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumBndNodes   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  NumEdges      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);  
  fclose(fids);
 
  % Read in the necessary "Cartesian" data (there will be duplicates here):
  fids  = fopen([outputdir,'/qhelp_cart.dat'],'r');
  if fids==-1
      error(['File  ',outputdir,'/qhelp_cart.dat  not found.']);
  end
  % assert that duplicates are the same:
  assert( meqn  == fscanf(fids,'%d',1) );
  assert( maux  == fscanf(fids,'%d',1) );
  assert( nplot == fscanf(fids,'%d',1) );
  assert( meth1 == fscanf(fids,'%d',1) );
  v_mx    = fscanf(fids,'%d',1);
  v_my    = fscanf(fids,'%d',1);
  v_xlow  = fscanf(fids,'%e',1);
  v_xhigh = fscanf(fids,'%e',1);
  v_ylow  = fscanf(fids,'%e',1);
  v_yhigh = fscanf(fids,'%e',1);
  assert( datafmt == fscanf(fids,'%e',1) );
  dvx = (v_xhigh-v_xlow)/v_mx;
  fclose(fids);

  % Number of polynomials in each direction:
  kmax2d = get_kmax(meth1, 2);
  kmax4d = get_kmax(meth1, 4);
  
%% ------------------------------------------------------------------------- %%
%% Create the mesh
%% ------------------------------------------------------------------------- %%

  disp(' Creating mesh ... ');
  
  %% READ-IN MESH INFO
  meshdir = [outputdir '/mesh_output'];
    
  fids = fopen([meshdir '/mesh_tnode.dat'],'r');
  tmp = fscanf(fids,'%d %d %d',[3,inf]);
  fclose(fids);
  tnode = transpose(tmp);
  clear tmp;
  tnode = tnode(1:NumPhysElems,1:3);

  fids = fopen([meshdir '/mesh_node.dat'],'r');
  tmp = fscanf(fids,'%e',[2,inf]);
  fclose(fids);
  node = transpose(tmp);
  clear tmp;
  node = node(1:NumPhysNodes,1:2);

%% ------------------------------------------------------------------------- %%
%% Default sampling data
%% ------------------------------------------------------------------------- %%

  xlow  = min(node(:,1));
  ylow  = min(node(:,2));
  xhigh = max(node(:,1));
  yhigh = max(node(:,2));
  
  xeps = 0.01*(xhigh-xlow);
  yeps = 0.01*(yhigh-ylow);   

%% ------------------------------------------------------------------------- %%
%% Select points for plotting:
%% ------------------------------------------------------------------------- %%

  % Find the information about the slice index of the grid:
  [BndyList, Psort, Mid, Left, Right] = FindSliceIndexUnst2( );

  % Add extra points and elements if points_per_dir>1
  node_old  = node;
  tnode_old = tnode;
% if (points_per_dir>1)
%   % This call modifies node (uncluding number of elements, etc.)!
%   [node, tnode, z2d_unst] = DivideUnst2Mesh(points_per_dir,  ...
%                                    NumPhysElems,node,tnode);    
% else
%   z4d_unst(1,1) = 0;
%   z4d_unst(1,2) = 0;
%   z4d_unst(1,3) = 0;
%   z4d_unst(1,4) = 0;
% end
  
  % Get physical midpoints of each element    
% xmid = zeros((NumPhysElems*points_per_dir^2),1);
% ymid = zeros((NumPhysElems*points_per_dir^2),1);
% for i=1:(NumPhysElems*points_per_dir^2)
%   xmid(i) = (node(tnode(i,1),1)+node(tnode(i,2),1)+ ...
%              node(tnode(i,3),1))/3;
%   ymid(i) = (node(tnode(i,1),2)+node(tnode(i,2),2)+ ...
%              node(tnode(i,3),2))/3;
% end
  
  % Sample Legendre polynomial on the midpoint of each element
% LegVals = GetHybr4Legendre(kmax4d, z4d_unst);
 
  disp(' Finished creating mesh. ');
  disp(' ');

  % Structure containing grid information (convienent for passing information):
  grid_struct.mx = v_mx;  
  grid_struct.my = v_my;
  grid_struct.NumElems     = length( BndyList );
  grid_struct.NumPhysElems = grid_struct.NumElems;
  %grid_struct.NumElems     = NumElems;
  %grid_struct.NumPhysElems = NumPhysElems;

  %% ---------------------------------- %%
  %% Construct a slice of the solution
  %% ---------------------------------- %%

  m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
              ' ) ? ']);
  disp(' ')
  if isempty(m)
    m=1;
  end
 
  % Flag for quitting
  q = -1;

  kn = 0;
  
  n  = 0;
  nf = 0;
  n1 = -1;
  while(nf~=-1)
    nf  = input([ ' Plot which frame ( 0 - ',num2str(nplot),...
                  ' ) [type -1 or q to quit] ? ']);
    if isempty(nf)
      n1 = n1 + 1;
      nf = 0;
    else
      n1 = nf;
    end
    if n1> nplot
      disp(' ');
      disp(' End of plots ');
      disp(' ');
      n1 = nplot;
    end
    if (nf~=-1)

      %% Solution -- q
      [q_data, time] = read_state4_hybrid(outputdir, n1, 'qslice', grid_struct, meqn, kmax4d);

      % Sample the function at vx=vy=0:
      qsoln = zeros( grid_struct.NumElems, grid_struct.mx, kmax4d );
      for k=1:kmax4d
      for vx=1:grid_struct.mx
      for ix=1:grid_struct.NumElems
        qsoln(ix, vx, k) = q_data( vx, ceil(grid_struct.my/2), ix, 1, k );
      end
      end
      end
      clear q_data;

%     if (maux>0)
%       %% Aux variables -- aux
%       [a_data,time] = read_state2_unst(datafmt, outputdir, n1, 'a', ...
%                                        NumElems, NumPhysElems, maux, kmax, 1:maux);
%       aux = sample_state2_unst(a_data, meth1, kmax, LegVals);
%       clear a_data;
%     end
      
      % USER SUPPLIED FUNCTION: Plotting function      
      NumPhysElems_new = NumPhysElems*points_per_dir^2;
      plotq4_hybrid_slices;
    end
    
  end
  disp(' ');
end
