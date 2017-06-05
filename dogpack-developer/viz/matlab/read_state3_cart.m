%function [q,time]=read_state3_cart(datafmt, basefilename, varname, ...
%  mx, my, mz, meqn, kmax, components)
function [q,time]=read_state3_cart(datafmt, outputdir, num_frame, varname, ...
                                   mx, my, mz, meqn, kmax, components)

  basefilename = [outputdir '/' varname sprintf('%04d',num_frame)];

  if(datafmt==1) % ASCII
    [q,time]=read_state3_ASCII(basefilename, varname, ...
                               mx, my, mz, meqn, kmax, components);
  else
    error(['unsupported value of datafmt: ' num2str(datafmt)]);
  end
  %toc;
end

function [q,time]=read_state3_ASCII(basefilename, varname, ...
                                    mx, my, mz, meqn, kmax, components)
  %
  filename = [basefilename '.dat'];
  fids = fopen(filename,'r');
  assert(fids~=-1, ['could not open file ' filename]);
  time = fscanf(fids,'%e',1);
  qsoln = fscanf(fids,'%e',[inf]);
  fclose(fids);  
  % order of subscripts must match WriteStateASCII in Output.cpp
  % (but reversed because matlab unlike C uses reverse odometer subscripts)
  %dims = [kmax,meqn,my,mx];
  dims = [mx,my,mz,meqn,kmax];

  q=reshape(qsoln,dims);
  % weed out the unused components
  if(~isequal(components,1:meqn))
    q=q(:,:,components,:);
  end
  % shift the subscripts to match HDF5 format
  permutation=[1,2,3,5,4];
  q=permute(q,permutation);
  %q=shiftdim(q,2);
end
