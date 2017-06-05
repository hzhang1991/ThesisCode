% plot field lines of divergenceless vector field
% by plotting contours of stream function
%
function plot_fieldlines(c, xc,yc,vecx,vecy,dx,dy,convect_fieldlines)
  delta_x = xc(2,1)-xc(1,1);
  delta_y = yc(1,2)-yc(1,1);
  assert(abs(dx-delta_x)<1e-15);
  assert(abs(dy-delta_y)<1e-15);
  stream_func = stream_func_from_vec(c, vecx,vecy,delta_x,delta_y);
  %stream_func = stream_func-stream_func(1,1);
  %min(min(yc))
  %max(max(yc))
  plot_contour(c, xc,yc,stream_func,convect_fieldlines);
end

function plot_contour(c, xc,yc,phi,convect_fieldlines)
   enforced_symmetry = c.s.params.enforced_symmetry;
   %phi=-phi;
   [mx,my]=size(xc);
   if(bitand(enforced_symmetry,3)==3)
      cx = [-xc(end:-1:1,end:-1:1),-xc(end:-1:1,:);...
            xc(:,end:-1:1),xc];
      cy = [-yc(end:-1:1,end:-1:1), yc(end:-1:1,:);...
            -yc(:,end:-1:1),yc];
      full_phi = [phi(end:-1:1,end:-1:1),phi(end:-1:1,:);...
                  phi(:,end:-1:1),phi];
      x_center = mx+1;
      y_center = my+1;
   elseif(bitand(enforced_symmetry,3)==2)
      cx = [ xc(:,end:-1:1), xc];
      cy = [-yc(:,end:-1:1), yc];
      full_phi = [phi(:,end:-1:1),phi];
   elseif(bitand(enforced_symmetry,3)==1)
      cx = [-xc(end:-1:1,:); xc];
      cy = [ yc(end:-1:1,:); yc];
      % assumes reflectional symmetry
      %full_phi = [phi(end:-1:1,:); phi];
      % assumes rotational symmetry
      full_phi = [phi(end:-1:1,end:-1:1); phi];
   else
      cx = xc;
      cy = yc;
      full_phi = phi;
   end
   plot_contour_helper(cx,cy,full_phi,convect_fieldlines);
end

function plot_contour_helper(cx,cy,full_phi,convect_fieldlines)
   assert(all(mod(size(full_phi),2)==[0 0]));
   %cy=cy-1.3;
   x_center = size(full_phi,1)/2+1;
   y_center = size(full_phi,2)/2+1;

   %maxp = max(max(phi));
   %minp = min(min(phi));
   %contour_step=(maxp-minp)/31;
   %upper=0:contour_step:maxp;
   %lower=0:-contour_step:minp;
   %levels=[lower(end:-1:2),upper];
   %levels=-15:15/15.;
   %contour(cx,cy,full_phi,levels,'-k');
   % should use a higher-order estimate of phi(0,0)
   %
   %pcolor(full_phi);
   hold on;
   num_field_lines = 24; % 16
   phi_center = full_phi(x_center,y_center); %+.000005
   %phi_center = interp2(full_phi,my+.5,mx+.5,'cubic');
   %contour(cx,cy,full_phi-phi_center-2,num_field_lines,'-k','linewidth',1);
   shifted_phi = full_phi-phi_center;
   if(convect_fieldlines)
     phi_min = full_phi(end,y_center); % this is only first-order accurate
     phi_max = full_phi(x_center,end); % this is basically second-order accurate
     %phi_max=max(max(full_phi));
     %phi_min=min(min(full_phi));
     phi_range=phi_max-phi_min;
     disp(['phi_range = ' num2str(phi_range) ' (should be 5.6900)']);
     delta_phi = phi_center - phi_min;
     %delta_phi=phi_center-phi_min;
     percent_reconnected = (delta_phi-.2)/(phi_range-.2);
     disp(['percent_reconnected = ' num2str(100*percent_reconnected) '%']);
     disp(['delta_phi = ' num2str(delta_phi)]);
     shifted_phi_min=phi_min-phi_center;
     %disp(['shifted_phi_min = ' num2str(shifted_phi_min)]);
     %percent_reconnected = -phi_min/phi_range;
     %disp(['percent_reconnected = ' num2str(percent_reconnected)]);
     contour_level_diff = phi_range/num_field_lines;
     contour_levels = (-num_field_lines:num_field_lines)*contour_level_diff;
     contour(cx,cy,shifted_phi+0.5*contour_level_diff,contour_levels,'-k','linewidth',1);
   else
     contour(cx,cy,full_phi,num_field_lines,'-k','linewidth',1);
   end
   contour(cx,cy,shifted_phi,[0,1000],'--r','linewidth',1);
end

function phi = stream_func_from_vec(c, B1,B2,dx,dy)
  mx=size(B1,1);
  my=size(B1,2);

  enforced_symmetry = c.s.params.enforced_symmetry;
  if(bitand(enforced_symmetry,1)==1)
    % define mesh representing integral of B2 from
    % each point's left neighbor to itself
    %B2_cell_xint = dx*0.5*padarray((B2(1:end-1,:)+B2(2:end,:)),[1,0],0,'pre');
    B2_cell_xint = dx*0.5*(B2(1:end-1,:)+B2(2:end,:));
    %
    % define mesh representing integral of B1 from
    % each point's upper neighbor to itself
    %B1_cell_yint = dy*padarray(0.5*(B1(:,1:end-1)+B1(:,2:end)),[0,1],'post');
    B1_cell_yint = dy*0.5*(B1(:,1:end-1)+B1(:,2:end));

    % for the sake of numerical stability we define the stream function to
    % be a weighted average of the path integral from the upper left to the given
    % point over all possible paths that move down and to the right along
    % lines between neighboring cell centers.  Inductively we say that the
    % value in each cell is an average of the value of the value in the
    % cell to the left plus the integral up to the present cell
    % with the value in the cell above plus the integral up to the present cell.
    %
    phi=zeros(size(B1));
    % first fill in the top boundary
    for ix=2:mx
      phi(ix,my) = phi(ix-1,my) + B2_cell_xint(ix-1,my);
    end
    % now work down the center
    for iy=my-1:-1:1
      phi(1,iy) = phi(1,iy+1) + B1_cell_yint(1,iy);
    end
    % now work down and to the right
    for iy=my-1:-1:1
    for ix=2:mx
      val_from_left = phi(ix-1,iy)+B2_cell_xint(ix-1,iy);
      val_from_abov = phi(ix,iy+1)+B1_cell_yint(ix,iy);
      phi(ix,iy) = 0.5*(val_from_left+val_from_abov);
    end
    end
  elseif(bitand(enforced_symmetry,1)==0)
    B2_cell_xint = dx*0.5*(B2(1:end-1,:)+B2(2:end,:));
    B1_cell_yint = dy*0.5*(B1(:,1:end-1)+B1(:,2:end));
    phi=zeros(size(B1));
    x_left_center = mx/2;
    x_rght_center = mx/2+1;
    % first fill in the top boundary
    for ix=(x_rght_center+1):mx
      phi(ix,my) = phi(ix-1,my) + B2_cell_xint(ix-1,my);
    end
    for ix=x_rght_center-1:-1:1
      phi(ix,my) = phi(ix+1,my) - B2_cell_xint(ix,my);
    end
    % now work down the (right) center
    for iy=my-1:-1:1
      phi(x_rght_center,iy) = ...
         phi(x_rght_center,iy+1) + B1_cell_yint(x_rght_center,iy);
    end
    % now work down and to the right
    for iy=my-1:-1:1
    for ix=x_rght_center+1:mx
      val_from_left = phi(ix-1,iy)+B2_cell_xint(ix-1,iy);
      val_from_abov = phi(ix,iy+1)+B1_cell_yint(ix,iy);
      phi(ix,iy) = 0.5*(val_from_left+val_from_abov);
      phi(ix,iy) = val_from_left;
    end
    end
    % now work down and to the left
    for iy=my-1:-1:1
    for ix=x_rght_center-1:-1:1
      val_from_rght = phi(ix+1,iy)-B2_cell_xint(ix+1,iy);
      val_from_abov = phi(ix,iy+1)+B1_cell_yint(ix,iy);
      phi(ix,iy) = 0.5*(val_from_rght+val_from_abov);
    end
    end
  end
end

