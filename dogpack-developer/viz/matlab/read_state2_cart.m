function [q, time] = read_state2_cart(datafmt, outputdir, num_frame, varname, ...
                                    mx, my, meqn, kmax, components, basefilename)
%READ_STATE2_CART.    Read in the 2D Cartesian state variable.
%
% This function samples the state variable saved in basefilename.  The default
% filename is given by
%
%     Default basefilename = outputdir/ [varname] + [num_frame] . [data_fmt],
%
% where data_type = ".dat" for ASCII and data_fmt = ".h5" for HDF5.
%
% This function serves the purpose of reformatting the file into a multi-index
% array with appropriate lengths.
%
% Usage: 
%
% [q, time] = read_state2_cart(datafmt, basefilename, varname, ...
%  mx, my, meqn, kmax, components)
%
% Input parameters:
%
%     datafmt - data format for output files
%             = 1:   ASCII format
%             = 5:   HDF5 format
%
%     outputdir - location of output directory
%
%     num_frame - frame number.  This goes into the file name that will be
%     read.
%
%     varname - name of the variable to be read.
%
%     mx, my - number of points in the x and y-direction.  This is important
%     for creating the proper size and dimensions for the state variable.
%
%     meqn - number of equations for the state vector being read in.
%
%     kmax - number of basis functions (per cell) for defining the state
%     vector
%
%     components - desired components to be used.  (Usually = 1:meqn).
%
%     basefilename - name of file, if different from
%         Default: outputdir/ [varname] + [num_frame] . [data_fmt],
%     where data_type = ".dat" for ASCII and data_fmt = ".h5" for HDF5.
%
% Output:
%
%     q - state vector.  size(q) = [mx, my, kmax, 1:len(components)].  
%
%     NOTE: this is NOT the same order of indexing as what happens in DoGPack.
%     The reason for this is to make the reader compatible with what happened
%     with the HDF5 part of the code.
%
%     time - time saved for this frame.  It is always assumed that this value
%     is saved as the first line in the output folder.
%
% See also: read_state2_unst, read_state3_cart.
    
  if(nargin<10)
      basefilename = [outputdir '/' varname ...
                      sprintf('%04d', num_frame)];
  end

  if(datafmt==1) % ASCII
    [q,time]=read_state2_ASCII(basefilename, varname, ...
                               mx, my, meqn, kmax, components);
  elseif(datafmt==5) % HDF5
    [q,time]=read_state2_HDF5(basefilename, varname, ...
                              mx, my, meqn, kmax, components);
  else
    error(['unsupported value of datafmt: ' num2str(datafmt)]);
  end

end

function [q,time]=read_state2_ASCII(basefilename, varname, ...
                                    mx, my, meqn, kmax, components)
%READ_STATE2_ASCII    Load in Cartesian data from an ASCII file.

  % Open the file and read in all of its contents:
  filename = [basefilename '.dat'];
  fids = fopen(filename,'r');
  assert(fids~=-1, ['could not open file ' filename]);

  % First line is always time, and every other line is a combination of basis
  % moments as well as equations at each cell in the grid:
  time  = fscanf(fids,'%e',1);
  qsoln = fscanf(fids,'%e',[inf]);
  fclose(fids);

  % Order of subscripts must match WriteStateASCII in Output.cpp
  % (but reversed because matlab unlike C uses reverse odometer subscripts)
  dims = [mx, my, meqn, kmax];
  q = reshape(qsoln, dims);

  % weed out the unused components (if any)
  if(~isequal(components, 1:meqn) )
    q=q(:,:,components,:);
  end

  % shift the subscripts to match HDF5 format  (why? -DS)
  permutation=[1, 2, 4, 3];
  q = permute(q, permutation);

end

function [q,time]=read_state2_HDF5(basefilename, varname, ...
                                   mx, my, meqn, kmax, components)
%READ_STATE2_HDF5    Load in Cartesian data from an HDF5 file.

  % In matlab it seems faster just to read the whole output
  % than to work with low-level access to array slices.
  filename = [basefilename '.h5'];
  fileinfo=hdf5info(filename);
  expected_dims=[mx,my,kmax,meqn];
  dims=fileinfo.GroupHierarchy.Datasets.Dims;
  assert(isequal(dims,expected_dims));
  time = fileinfo.GroupHierarchy.Datasets.Attributes.Value;
  % it is apparently more efficient to read in the entire data;
  q=hdf5read(filename,varname);
  [d1,d2,d3,d4]=size(q);
  disp(['finished reading HDF5 array with dimensions ' num2str(size(q))]);
  assert(all(dims==[d1,d2,d3,d4]));
  if(~isequal(components,1:meqn))
    q=q(:,:,:,components);
  end

end
