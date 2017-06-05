
function initial_yaxis_flux = get_initial_yaxis_flux(xplotter,outputdir)

  % read first line of reconnected flux file
  %
  fid = fopen([outputdir '/recon_flux.dat']);
  consv = fscanf(fid,'%e',[4 1]);
  status = fclose(fid);

  initial_yaxis_flux = consv(4,1);
end

