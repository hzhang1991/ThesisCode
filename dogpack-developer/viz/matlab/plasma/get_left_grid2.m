
function [dx,dy,xl,yl]=get_left_grid2(mx,my,xlow,xhigh,ylow,yhigh,method1)
  dx = (xhigh-xlow)/mx;
  dy = (yhigh-ylow)/my;
  % low/left edges of cells
  xl=xlow+(0:mx*method1)*(dx/method1);
  yl=ylow+(0:my*method1)*(dy/method1);
  % should the user do this?
  [xl,yl]=ndgrid(xl,yl);
end

