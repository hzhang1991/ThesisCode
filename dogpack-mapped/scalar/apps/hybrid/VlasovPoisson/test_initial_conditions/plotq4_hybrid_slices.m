%% -------------------------------------------------------------------------- %%
%% Scatter plot of the exact Solution
%% -------------------------------------------------------------------------- %%

NumPhysElems_big = NumPhysElems;
NumPhysElems     = grid_struct.NumPhysElems;

%% -------------------------------------------------------------------------- %%
%% Construct a slice of the solution at y = 0
%% -------------------------------------------------------------------------- %%


% Grab the slice index for this problem:
[BndyList, Psort, MidPt, LftPt, RghtPt] = FindSliceIndexUnst2( );
NumTriOnLine = length( BndyList ); 


%   % Sample the function on each triangle (this won't be a uniform sampling ... )
%   [q_data, T] = read_state2_unst(datafmt, outputdir, n1, 'q', ...
%       NumElems, NumPhysElems, meqn, kmax4d, 1:meqn);      

local_pts_per_dir = points_per_dir;
rline = zeros( local_pts_per_dir*NumTriOnLine, 1 );
qline = zeros( local_pts_per_dir*NumTriOnLine, 1 );
index = 0;
for ii=1:NumTriOnLine

      % Find the coefficients for this triangle:
      TriNum = BndyList( Psort(ii) );
%     qtmp   = zeros(1, meqn, kmax4d);
%     qtmp(1, :, : ) = q_data( TriNum, : , : );

      % Center of current triangle:
      % REMEMBER TO USE OLD HERE: DivideUnst2Mesh modifies node and tnode!!!
      x_tmp = node_old( tnode_old(TriNum,:), 1 );   x_m = sum( x_tmp ) / 3;
      y_tmp = node_old( tnode_old(TriNum,:), 2 );   y_m = sum( y_tmp ) / 3;

      % Triangle boundary
      x1 = x_tmp(1);  y1 = y_tmp(1);
      x2 = x_tmp(2);  y2 = y_tmp(2);
      x3 = x_tmp(3);  y3 = y_tmp(3);
      A = [ (x2-x1), (x3-x1); (y2-y1), (y3-y1) ];

%     % Canonical points we wish to sample:
      s2d_tmp = zeros( local_pts_per_dir, 2 );
      x = zeros( local_pts_per_dir, 1 );
      xl_tmp = LftPt( Psort(ii), 1 );  xr_tmp = RghtPt( Psort(ii), 1 );
      for mi=1:local_pts_per_dir
        x(mi)         = xl_tmp + (mi-0.5)*( xr_tmp-xl_tmp )/local_pts_per_dir;
%       assert( x(mi) >= xl_tmp & x(mi) <= xr_tmp );
        s2d_tmp(mi,:) = A \ [ x(mi)-x_m; -y_m ];
      end
      s4d_tmp = zeros( local_pts_per_dir, 4 );
      s4d_tmp(:,1:2) = s2d_tmp(:,1:2);

      % Need to reconstruct LegVals each time (because it will be different
      % for each cell).
      LegVals_tmp = GetHybr4Legendre(kmax4d, s4d_tmp);
      qvals       = qsoln( Psort(ii), : ) * LegVals_tmp;

      for mp=1:local_pts_per_dir
        index = index+1;
        rline(index) = x(mp);
        qline(index) = qvals(mp);
      end
end
clear q_data;

%% -------------------------------------------------------------------------- %%
%% Plots
%% -------------------------------------------------------------------------- %%

% Single plot of a slice through y = 0
figure(1);
clf;
pz = plot( rline, qline, 'bo' );
set(pz, 'linewidth', 2);
axis on; box on; grid off;
axis([-1 1.0 -0.1 1.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(t,y=0) at t = ',...
    num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
hold on;

% Include plot of exact solution:
%p2 = plot( rline, (0.75+0.25*cos(2*pi*time)) * (0.5+0.5*cos(pi*rline)), 'r-' ); 

figure(1)

NumPhysElems = NumPhysElems_big;
