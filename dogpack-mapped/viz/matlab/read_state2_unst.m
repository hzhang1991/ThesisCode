function [q, time] = read_state2_unst(datafmt, outputdir, num_frame, varname, ...
                                   NumElems, NumPhysElems, meqn, kmax, components)
%READ_STATE2_UNST.    Read in a 2D unstructured state variable
%
% This function reads in a variable labeled by varname from the output
% directory ourputdir, and returns a multi-dimensional array allowing easy
% indexing into that variable.
%
% Usage: 
%
% [q, time] = read_state2_unst(datafmt, basefilename, varname, ...
%                                      mx, my, meqn, kmax, components)
%
% The format of the file to be read in is as follows:
%
%     basefilename = outputdir/ [varname] + [num_frame] . [data_fmt],
%     where data_type = ".dat" for ASCII and data_fmt = ".h5" for HDF5.
%     (Currently only supports ASCII type data files).
%
% Input parameters:
%
%     datafmt - data format for output files  (currently unused!)
%             = 1:   ASCII format  :  This is currently the only format that
%             is read
%             = 5:   HDF5 format   : Future versions will allow for this.
%
%     outputdir - location of output directory
%
%     num_frame - frame number.  This goes into the file name that will be
%     read.
%
%     varname - name of the variable to be read.
%
%     NumElems, NumPhysElems - Grid information.  This is important
%     for creating the proper size and dimensions for the state variable.
%     These parameters need to identically match the grid information used to
%     produce the solution.
%
%     meqn - number of equations for the state vector being read in.
%
%     kmax - number of basis functions (per cell) for defining the state
%     vector
%
%     components - desired components to be saved.  (Usually = 1:meqn).  In
%     the case of a solution with many equations (e.g. high-moment two-fluid)
%     one may only want to sample and plot a few of these equations.
%
% Output:
%
%     q - state vector.  size(q) = [NumPhysElems, 1:len(components), kmax].
%
%     time - time saved for this frame.  It is always assumed that this value
%     is saved as the first line in the output folder.
%
% See also: read_state2_cart, read_state3_cart

  basefilename = [outputdir '/' varname sprintf('%04d',num_frame)];
  filename = [basefilename '.dat'];
  fids  = fopen(filename,'r');
  assert(fids~=-1, ['could not open file ' filename]);
  time  = fscanf(fids,'%e',1);
  qsoln = fscanf(fids,'%e',[inf]);
  fclose(fids);  

  dims = [NumElems,meqn,kmax];
  q    = reshape(qsoln,dims);
  
  % weed out the unused components
  q = q(1:NumPhysElems,components,:);

end
