
% could implement this with a call to read_config_file
%
function [gamma,B_0,BCs,domain_scaling,enforced_symmetry, ...
  mass_ratio,temp_ratio,cc,cs_light,clean_E_field] ...
  = read_param_plasma(outputdir)

  readfile=[outputdir '/param.data'];

  global gamma B_0 BCs domain_scaling enforced_symmetry ...
    mass_ratio temp_ratio cc cs_light clean_E_field;
  % set default values
  domain_scaling=1;
  enforced_symmetry=7; %(1=X)+(2=Y)+(4=actual_mesh)
  % read the values
  read_config_file(readfile);
end
