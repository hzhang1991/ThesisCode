function CreateMatrix(morder,showplots)
%CREATEMATRIX   Create the matrix for Poisson's equation.
%
% This wrapper function calls CreateMatrixGlobal.  This function forces
% periodic boudary conditions to be defined.
%
% Usage: CREATEMATRIX( morder, showplots )
%
% Inputs:
%
% morder    = 1 (CG1 method)
%           = 2 (CG2 method -- needs to match SubFactor parameter)
%           = 3 (CG2 method -- needs to match SubFactor parameter)
%
% showplots = 1 (default: show sparsity plots of matrices)
%             0 (don't show sparsity plots)
%
% See also: CreateMatrixGlobal

%   addpath('../../lib/matlab');
    bctype=2;
    CreateMatrixGlobal(morder,bctype,showplots);

end
