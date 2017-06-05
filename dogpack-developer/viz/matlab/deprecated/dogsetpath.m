
% ensure that MATLABPATH includes the expected 
% directories in the expected order
%
  DOGPACK = getenv('DOGPACK');
  assert(~isempty(DOGPACK),...
    ['Must set DOGPACK environment variable to point to DOGPACK directory']);

  addpath(...
     [DOGPACK '/matlab'],...
     [DOGPACK '/matlab/deprecated']);

