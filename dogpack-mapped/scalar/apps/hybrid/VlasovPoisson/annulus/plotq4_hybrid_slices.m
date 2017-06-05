%% -------------------------------------------------------------------------- %%
%% Scatter plot of the exact Solution
%% -------------------------------------------------------------------------- %%

NumPhysElems_big = NumPhysElems;
NumPhysElems     = grid_struct.NumPhysElems;

mx = grid_struct.mx;

%% -------------------------------------------------------------------------- %%
%% Construct a slice of the solution at y = 0
%% -------------------------------------------------------------------------- %%

% DO NO REMOVE THIS ASSERT STATEMENT UNLESS YOU KNOW WHAT YOU ARE DOING.  This
% will break everything below if removed.
assert( points_per_dir == 1 );

% Grab the slice index for this problem:
[BndyList, Psort, MidPt, LftPt, RghtPt] = FindSliceIndexUnst2( );
NumTriOnLine = length( BndyList ); 

%   % Sample the function on each triangle (this won't be a uniform sampling ... )
%   [q_data, T] = read_state2_unst(datafmt, outputdir, n1, 'q', ...
%       NumElems, NumPhysElems, meqn, kmax4d, 1:meqn);      

local_pts_per_dir = points_per_dir;
npts = local_pts_per_dir^2;

xx    = zeros( local_pts_per_dir*NumTriOnLine, local_pts_per_dir*mx );
yy    = zeros( local_pts_per_dir*NumTriOnLine, local_pts_per_dir*mx );
qpts  = zeros( local_pts_per_dir*NumTriOnLine, local_pts_per_dir*mx );

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

      % Canonical points we wish to sample:
      s2d_tmp = zeros( local_pts_per_dir, 2 );
      x = zeros( local_pts_per_dir, 1 );
      xl_tmp = LftPt( Psort(ii), 1 );  xr_tmp = RghtPt( Psort(ii), 1 );
      for mi=1:local_pts_per_dir
        x(mi)         = xl_tmp + (mi-0.5)*( xr_tmp-xl_tmp )/local_pts_per_dir;
        s2d_tmp(mi,:) = A \ [ x(mi)-x_m; -y_m ];
      end
      s4d_tmp = zeros( local_pts_per_dir^2, 4 );
      s4d_tmp(:,1:2) = s2d_tmp(:,1:2);

      % UNIFORM points in velocity space!  Common to every element
      s4d_tmp(:,  3) = -1+(2*(1:points_per_dir)-1)/points_per_dir;

      % Need to reconstruct LegVals each time (because it will be different
      % for each triangular cell).  We can recycle these for each element in
      % Cartesian space.
      LegVals_tmp = GetHybr4Legendre(kmax4d, s4d_tmp);

      tmp_soln = zeros( 1, kmax4d );
      for vx=1:mx
        % Sample at the current point:
        tmp_soln(:) = qsoln( Psort(ii), vx, 1:kmax4d );
        qvals  = tmp_soln * LegVals_tmp;
        for mp=1:local_pts_per_dir
          xx(index+mp, vx) = x(mp);
          yy(index+mp, vx) = v_xlow + (vx-0.5)*dvx;
          qpts  (index+mp, vx) = qvals(mp);
        end
      end
      index = index+local_pts_per_dir;

end
clear q_data;

%% -------------------------------------------------------------------------- %%
%% Plots
%% -------------------------------------------------------------------------- %%

% TODO - use these units!
% units with meaning:
N = 1.0e13; % 1/m
L = 0.1; % m
T = 5.605424055e-9; % s
E0 = 1.809512620e4; % N/C
Phi0 = 1.809512620e3; % Nm/C
%  F = 5.605424055e5; % s/m^2

% Single plot of a slice through y = 0
qmx_axis = 3.0e15;
figure(1);
clf;
pz = plot( L*xx(:, ceil(mx/2)), N*qpts(:, ceil(mx/2)), 'bo' );
set(pz, 'linewidth', 2);
axis on; box on; grid off;
axis([-0.1 0.1 -0.1 qmx_axis]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-0.1:0.02:0.1);
set(gca,'ytick',0:0.25e15:qmx_axis);
set(gca,'fontsize',16);
t1 = title(['q(t,x,y=0,vx= ', num2str( yy(1,ceil(mx/2))), ' vy=0) at t = ',...
       num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
hold on;

figure(2)
clf;

  pcolor(xx, yy, qpts )

  colormap('jet');
  axis on; box on; grid off;
  set(gca,'fontsize',16);
  t1 = title(['f(t,x,v) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);
  shading flat;
  c1 = colorbar;
  caxis auto;
% set(gca,'fontsize',16);
% set(gca,'xtick',-6:2:6);
% set(gca,'ytick',-6:2:6);


%   NumPhysElems = NumPhysElems_big;
