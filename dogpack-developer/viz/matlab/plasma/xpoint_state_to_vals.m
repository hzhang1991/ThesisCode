% for second-order data all this really does is to extract
% the zeroth coefficient
%
function xpoint_state_to_vals(outputdir_in)
  format long e
  setarg_outputdir;
  parameters_ini = [outputdir '/parameters.ini'];
  params = read_params_from_ini(parameters_ini, 'meqn', 'space_order');
  meqn = params.meqn;
  kmax = get_kmax(params.space_order);
  
  fid = fopen([outputdir '/xpoint_cell_state.dat']);
  xpoint_cell_state = fscanf(fid,'%e',[meqn*kmax+1, inf]);
  numTimeSteps = size(xpoint_cell_state,2);
  fclose(fid);
  q = reshape(xpoint_cell_state(2:end,:), [kmax, meqn, numTimeSteps]);
  phi = sample_basis_functions2(params.space_order, 1, 1);
  assert(all(size(phi)==[1,1,kmax]));
  qvals = zeros(1,meqn, numTimeSteps);
  for k=1:kmax
    qvals=qvals+phi(1,1,k)*q(k,:,:);
  end
  fid=fopen([outputdir '/xpoint_cell_vals.dat'],'w');
  for timeStep=1:numTimeSteps
    fprintf(fid,'%.16e ', xpoint_cell_state(1,timeStep));
    fprintf(fid,'%.16e ', qvals(1,:,timeStep));
    fprintf(fid,'\n');
  end
  fclose(fid);

end
