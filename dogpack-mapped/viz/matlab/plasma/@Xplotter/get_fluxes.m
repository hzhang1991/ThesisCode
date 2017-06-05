function [fluxes] = get_fluxes(xplotter, outputdir)
  
  %suffix = '.16.dat';
  suffix = '.dat';
  
  fid = fopen([outputdir '/recon_flux' suffix]);
  consv = fscanf(fid,'%e',[6 inf]);
  status = fclose(fid);

  initial_yaxis_flux = consv(4,1);

  fluxes.t = consv(1,:)';
  fluxes.change_in_bottom_to_rght_flux = consv(2,:)'-consv(2,1);
  fluxes.rght_flux = consv(3,:)';
  fluxes.left_flux = consv(4,:)';
  fluxes.recon_flux = fluxes.left_flux(1,1) - fluxes.left_flux;
  fluxes.left_lost_flux = consv(5,:)';
  total_x_axis_island_flux = consv(6,:)';

  fluxes.show_central_island_flux = 0;
  if(max(abs(fluxes.left_lost_flux))/initial_yaxis_flux > 0.01)
    fluxes.show_central_island_flux = 1;
  end
  %
  fluxes.other_island_flux = total_x_axis_island_flux-fluxes.left_lost_flux;
  fluxes.show_other_island_flux = 0;
  if(max(abs(fluxes.other_island_flux))/initial_yaxis_flux > 0.01)
    fluxes.show_other_island_flux = 1;
  end
  
  fid = fopen([outputdir '/xpoint_E3' suffix]);
  if(fid ~= -1)
    xpoint_E3 = fscanf(fid,'%e',[2 inf]);
    status = fclose(fid);
    %
    E3 = xpoint_E3(2,:)';
    E3t = xpoint_E3(1,:)';
    %if(length(E3t)~=length(fluxes.t));
    %  warning(['mismatch: length(E3t)=' num2str(length(E3t)) ...
    %    ' but length(fluxes.t)=' num2str(length(fluxes.t))]);
    %  len=min([length(E3t) length(fluxes.t)]);
    %  E3t = E3t(1:len);
    %  t = fluxes.t(1:len);
    %  indices = find(E3t~=t);
    %  if(numel(indices)>0)
    %    idx = indices(1);
    %    idx_s = num2str(idx);
    %    warning(['mismatch: E3t(' idx_s ') = ' num2str(E3t(idx)) ...
    %      ' but t(' idx_s ') = ' num2str(t(idx))]);
    %  end
    %end
    assert(all(fluxes.t==E3t));
    assert(length(E3)==length(fluxes.t));
    fluxes.intE3 = integrate(E3t,E3);
  end
end

function out=integrate(t,f)
  % calculute an accumulation integral using
  % the trapezoid rule
  fh = (f(2:end) + f(1:end-1))*0.5;
  dt = t(2:end)-t(1:end-1);
  out=zeros(size(f));
  out(2:end) = cumsum(fh.*dt);
end
