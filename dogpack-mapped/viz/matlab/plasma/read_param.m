
% could implement this with a call to read_config_file
%
function [gamma,B_0,BCs,domain_scaling,cc,enforced_symmetry] = read_param(outputdir)

  readfile=[outputdir '/param.data'];

  global gamma B_0 BCs domain_scaling cc enforced_symmetry;
  % set default values
  domain_scaling=1;
  enforced_symmetry=3;
  % read the values
  read_config_file(readfile);
end
