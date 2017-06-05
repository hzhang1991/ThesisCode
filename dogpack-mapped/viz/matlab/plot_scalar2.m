
function plot_scalar2(fig, values, xl, yl, time, axis_array, name);
  % must augment arrays to accomodate ideosyncrasies of
  % matlab plotting routines.
  values_aug = zeros(size(values)+1);
  values_aug(1:end-1,1:end-1) = values;
  xl_full = kron(xl',ones(1,numel(yl)));
  yl_full = kron(ones(numel(xl),1),yl);
  figure(fig);
  clf;
  % pcolor and surf seem to give identical behavior,
  % except that surf leaves a faint white border on the top and right..
  pcolor(xl_full,yl_full,values_aug);
  %pcolor(xl_full(1:end-1,1:end-1),yl_full(1:end-1,1:end-1),values_aug(1:end-1,1:end-1));
  axis('equal');
  axis(axis_array);
  set(gca,'fontsize',16);
  t1 = title([name ' at t = ', num2str(time), '     [pvFrame]']); 
  set(t1,'fontsize',16);
  c1=colorbar;
  set(c1,'fontsize',16);
  shading flat;
end
