
function [t, rho_i, M3i, P13i_x, P23i_y, E3] = read_xpoint_data(xplotter,...
  outputdir, species)
  %
  suffix = '.dat';
  %
  readfile = [outputdir '/xpoint_gas_' species suffix];
  disp(['reading from file ' readfile]);
  fid = fopen(readfile);
  num_columns = numel(sscanf(fgetl(fid),'%e',[inf]));
  assert(num_columns>=5);
  frewind(fid);
  xpoint_gas = fscanf(fid,'%e',[num_columns inf]);
  status = fclose(fid);
  %
  t = xpoint_gas(1,:)';
  rho_i = xpoint_gas(2,:)';
  M3i = xpoint_gas(3,:)';
  P13i_x = xpoint_gas(4,:)';
  P23i_y = xpoint_gas(5,:)';
  
  fid = fopen([outputdir '/xpoint_E3' suffix]);
  xpoint_E3 = fscanf(fid,'%e',[2 inf]);
  status = fclose(fid);
  %
  E3 = xpoint_E3(2,:)';
  if(length(E3)~=length(t));
    warning(['mismatch: length(E3)=' num2str(length(E3)) ...
      ' but length(t)=' num2str(length(t))]);
  end
end

