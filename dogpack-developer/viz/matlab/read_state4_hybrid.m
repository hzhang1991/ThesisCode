function [q,time]=read_state4_hybrid(outputdir,num_frame,varname,grid,meqn,kmax)
%READ_STATE4_HYBRID.     Read from file the state vector.
%
% Usage:
%
%  [q,time] = READ_STATE4_HYBRID( outputdir,num_frame,varname,grid,meqn,kmax )
%
% Reads in the state vector from file and reshapes its size for easier
% indexing.  The state vector in file is a single long list of digits that
% makes for difficult indexing.  This routine turns it into a multi-d array
% for easier indexing.
%
% Parameters:
%
%    outputdir - location of the output directory
%
%    varname, num_frame - used to identify what file to read.
%
%    grid - struct containing (at least) the following elements:
%           { mx, my, NumElems, NumPhysElems }
%
%    meqn - number of equations
%    kmax - number of polynomials
%
% Output:
%
%    q    - state vector.  size(q) = [mx,my,NumElems,meqn,kmax]
%
%    time - time of solution that is read.  
%           Current format for files is this
%           is saved in the first line of the output file.
%
% See also: SAMPLE_STATE4_HYBRID

    mx = grid.mx; my = grid.my;
    NumElems     = grid.NumElems;
    NumPhysElems = grid.NumPhysElems;

    % Scan in the filename.  For the order in which elements are saved, see, for
    % example: $(DOGPACK)/scalar/lib/Ouput_Hybrid.cpp.
    %
    % The ouput format is a single long list of numbers, the first
    % element is time.
    basefilename = [outputdir '/' varname sprintf('%04d',num_frame)];
    filename = [basefilename '.dat'];
    fids = fopen(filename,'r');
    assert(fids~=-1, ['could not open file ' filename]);
    time  = fscanf(fids,'%e',1);
    qsoln = fscanf(fids,'%e',[inf]);
    fclose(fids);  

    % Reshape the solution into a multi-D array with indices that have meaning:
    dims = [mx, my, NumElems, meqn, kmax];
    q    = reshape(qsoln, dims);

end
