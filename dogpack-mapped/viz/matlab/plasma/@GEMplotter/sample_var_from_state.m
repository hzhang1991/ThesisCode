
% For documentation on how to implement these sampling functions,
% see sample_dtJ for out-of-date documentation on data representation.
% Perhaps a better way to handle this would be:
% + create a "class" (cell array) which contains
%   - a data format indicator field and
%   - the array
% + convert format as needed
% + pass to each sample_* function only the data it requires
% + unix philosophy: split up into small routines each of which
%   does one thing well.
%
function var_vals  = sample_var_from_state(c, var_name, state, varargin)

  % set default keyword arguments
  do_rescale=1;
  %
  % parse keyword arguments
  %
  idx=1;
  while idx <= numel(varargin)
    switch lower(varargin{idx})
      case 'do_rescale'; idx=idx+1;
        do_rescale=varargin{idx}; idx=idx+1;
      otherwise
        error(['invalid keyword: ' varargin{idx}]);
    end
  end
  
  if(isfield(c.s.stateVarMap.varIdx_s,var_name))
    var_vals  = sample_state_var_from_state(c, var_name, state);
  elseif(isfield(c.s.varMap.varIdx_s,var_name))
  % This is a hack; I should really use an array of function pointers
  % rather than generate code.
    command = ['var_vals = sample_' var_name '(c, state);'];
    eval(command);
  else
    error(['no such variable: ' var_name]);
  end

  % rescale the data values if requested
  %
  var_scale=1;
  if(do_rescale)
    varIdx = c.s.varMap.varIdx_s.(var_name);
    var_scale = c.s.plot_params.output_scales(:,varIdx);
    var_scale;
  end
  %var_scale
  if(var_scale~=1)
    disp(['dividing ' var_name ' values by factor of ' num2str(var_scale)]);
    var_vals = var_vals./var_scale;
  end
end

function dims = get_plot_arr_dims(c)
  if(c.s.sample_rate==0)
    numpts = c.s.params.space_order;
  else
    numpts = c.s.sample_rate;
  end
  plot_mx=c.s.params.plot_mx;
  plot_my=c.s.params.plot_my;
  mx_out = numpts*plot_mx;
  my_out = numpts*plot_my;
  dims = [mx_out,my_out];
end

function arr = sample_zeros(c,meqn)
  dims = get_plot_arr_dims(c);
  arr=zeros([dims,meqn]);
end

function arr = sample_ones(c,meqn)
  dims = get_plot_arr_dims(c);
  arr=ones([dims,meqn]);
end

function var_vals  = sample_state_var_from_state(c, var_name, state_coef, sample_rate)
    if(~exist('sample_rate','var')); sample_rate = c.s.sample_rate; end;
    components = c.s.stateVarMap.stateIndices_s.(var_name);
    var_coef = state_coef(:,:,:,components);
    var_vals = sample_state_cart2(var_coef, sample_rate);
end

%%%%%%%%%%%%%%%%%%%%%%
%%% basic routines %%%
%%%%%%%%%%%%%%%%%%%%%%

function out = val_crossProduct(A,B)
  assert(all(size(A)==size(B)));
  assert(size(size(A),2)==3);
  assert(any(size(A,3)==[2 3]));
  out = zeros(size(A));
  out(:,:,1) = A(:,:,2).*B(:,:,3) - A(:,:,3).*B(:,:,2);
  out(:,:,2) = A(:,:,3).*B(:,:,1) - A(:,:,1).*B(:,:,3);
  out(:,:,3) = A(:,:,1).*B(:,:,2) - A(:,:,2).*B(:,:,1);
end

% compute the third component of the cross product
%
function out = val_crossProduct3(A,B)
  assert(all(size(A)==size(B)));
  assert(size(size(A),2)==3);
  assert(any(size(A,3)==[2 3]));
  out = zeros(size(A,1),size(A,2),1);
  out(:,:,1) = A(:,:,1).*B(:,:,2) - A(:,:,2).*B(:,:,1);
end

function product = tupleOverScalar(tuple,scalar)
  product = zeros(size(tuple));
  for idx=1:size(tuple,3)
    product(:,:,idx) = tuple(:,:,idx)./scalar;
  end
end

function product = tupleTimesScalar(tuple,scalar)
  product = zeros(size(tuple));
  for idx=1:size(tuple,3)
    product(:,:,idx) = tuple(:,:,idx).*scalar;
  end
end


% get eigenvalues of symmetric tensor
% (eigenvalues appear to be in order from least to greatest)
function [eig1, eig2, eig3] = get_eigs(P11, P12, P13, P22, P23, P33);
  %
  a1=-P11-P22-P33;
  a2=P11.*P22+P11.*P33+P22.*P33-P12.^2-P13.^2-P23.^2;
  a3=P11.*P23.^2+P22.*P13.^2+P33.*P12.^2 ...
    -2*P12.*P13.*P23-P11.*P22.*P33;
  Q=(3*a2-a1.^2)/9;
  R=(9*a1.*a2-27*a3-2*a1.^3)/54;
  discriminant=Q.^3+R.^2;
  % discriminant should be negative
  %find(discriminant>0)
  root_disc = sqrt(discriminant); % should be imaginary
  S=(R+root_disc).^(1/3);
  % could just let T be the complex conjugate of S.
  T=(R-root_disc).^(1/3);
  % should be real
  S_plus_T=real(S+T);
  imag_S_minus_T = imag(S-T);
  % The vast majority of the work is necessary to find even
  % one eigenvalue.
  eig3=S_plus_T-a1/3;
  partA = -S_plus_T/2-a1/3;
  partB = (sqrt(3)/2)*imag_S_minus_T;
  eig2= partA+partB;
  eig1= partA-partB;
end

% accessing state variables

function m_s = get_species_mass(c, species)
  switch species
  case 'i'
    m_s = c.s.params.ion_mass;
  case 'e'
    m_s = c.s.params.elc_mass;
  end
end

function model=get_species_model(c, species);
 model = 5;
 switch c.s.params.model_name
 case 'i10e5'
  if(strcmp(species,'i')) model=10; end
 case 'g10'
  model=10;
 case 'p10'
  model=10;
 end
end

function [rho, u1, u2, u3, p] = get_primVars05(c, state_coef, species);
  model=get_species_model(c, species);

 if(model==10)
  [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
     = get_primVars10(c, state_coef, species);
  p = (P11 + P22 + P33)/3.;
 else % model==5
  [rho_s_, Ms_, nrg_s_] = get_gasVarNames05(species);

  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  % get state coefficients
  rho_s_coef = state_coef(:,:,:,stateIndices_s.(rho_s_));
  Ms1_coef = state_coef(:,:,:,stateIndices_s.(Ms_)(1));
  Ms2_coef = state_coef(:,:,:,stateIndices_s.(Ms_)(2));
  Ms3_coef = state_coef(:,:,:,stateIndices_s.(Ms_)(3));
  nrg_s_coef = state_coef(:,:,:,stateIndices_s.(nrg_s_));
  % get kinetic energy
  rho  = sample_state_cart2(rho_s_coef,c.s.sample_rate);
  M1 = sample_state_cart2(Ms1_coef,c.s.sample_rate);
  M2 = sample_state_cart2(Ms2_coef,c.s.sample_rate);
  M3 = sample_state_cart2(Ms3_coef,c.s.sample_rate);
  u1 = M1./rho;
  u2 = M2./rho;
  u3 = M3./rho;
  nrg = sample_state_cart2(nrg_s_coef,c.s.sample_rate);
  KE = (u1.*M1 + u2.*M2 + u3.*M3)/2.;
  % get pressure
  p = (c.s.params.gamma-1.)*(nrg-KE);
 end
end

function [eig1, eig2, eig3, rho, u1, u2, u3] ...
  = sample_eigs_s(c, state_coef, species)

  species_model=get_species_model(c, species);
  switch species_model
    case 5
      [rho, u1, u2, u3, p] = get_primVars05(c, state_coef, species);
      eig1 = p;
      eig2 = p;
      eig3 = p;
    case 10
      [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
         = get_primVars10(c, state_coef, species);

      % compute eigenvalues of total pressure tensor
      [eig1, eig2, eig3] = get_eigs(P11, P12, P13, P22, P23, P33);
      % verify that the pressure tensor is positive definite.
      if(0==all(all(eig1>0.)))
          disp(['eig1 is not everywhere positive for species ' species]);
      end
      if(0==all(all(eig2>0.)))
          disp(['eig2 is not everywhere positive for species ' species]);
      end
      if(0==all(all(eig3>0.)))
          disp(['eig3 is not everywhere positive for species ' species]);
      end
      % assert(all(all(eig2>0.)));
      % assert(all(all(eig3>0.)));

    otherwise
      error(['unsupported value: ' num2str(species_model)]);
  end
end

function [rho_s, Ms, Ns] = get_gasVarNames10(species);
  rho_s = ['rho_' species];
  Ms    = ['M' species];
  Ns    = ['N' species];
end

function [rho_s, Ms, nrg_s] = get_gasVarNames05(species);
  rho_s = ['rho_' species];
  Ms    = ['M' species];
  nrg_s = ['nrg_' species];
end

function [rho_s, Ms] = get_gasVarNames4(species);
  rho_s = ['rho_' species];
  Ms    = ['M' species];
end

function [rho_s, Ms, energy_s] = get_gasVarNames(c, species);
  species_model=get_species_model(c, species);
  switch species_model
  case 10
    [rho_s, Ms, energy_s] = get_gasVarNames10(species);
  case 5
    [rho_s, Ms, energy_s] = get_gasVarNames05(species);
  otherwise
    error(['unsupported model: ' num2str(species_model)]);
  end
end

function [rho, M1, M2, M3] = get_consVars4(c, state_coef, species)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  [rho_s_, Ms_] = get_gasVarNames4(species);
  rho_coef = state_coef(:,:,:,stateIndices_s.(rho_s_));
  M1_coef  = state_coef(:,:,:,stateIndices_s.(Ms_)(1));
  M2_coef  = state_coef(:,:,:,stateIndices_s.(Ms_)(2));
  M3_coef  = state_coef(:,:,:,stateIndices_s.(Ms_)(3));
  rho = sample_state_cart2(rho_coef,c.s.sample_rate);
  M1 = sample_state_cart2(M1_coef,c.s.sample_rate);
  M2 = sample_state_cart2(M2_coef,c.s.sample_rate);
  M3 = sample_state_cart2(M3_coef,c.s.sample_rate);
end

function [rho, u1, u2, u3, M1, M2, M3] = get_primVars4(c, state_coef, species)
  [rho, M1, M2, M3] = get_consVars4(c, state_coef, species);
  rho_inv = 1./rho;
  u1 = M1.*rho_inv;
  u2 = M2.*rho_inv;
  u3 = M3.*rho_inv;
end

% species: 'i' or 'e'
function [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
  = get_primVars10(c, state_coef,species)

  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  %
  % get state coefficients
  %
  Ns_ = ['N' species];
  N11_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(1));
  N12_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(2));
  N13_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(3));
  N22_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(4));
  N23_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(5));
  N33_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(6));

  % get kinetic energy tensor
  %
  [rho, u1, u2, u3, M1, M2, M3] = get_primVars4(c, state_coef, species);
  K11 = M1.*u1;
  K12 = M1.*u2;
  K13 = M1.*u3;
  K22 = M2.*u2;
  K23 = M2.*u3;
  K33 = M3.*u3;
  %
  % get elc pressure state coefficients
  %
  % access energy tensor
  N11 = sample_state_cart2(N11_coef,c.s.sample_rate);
  N12 = sample_state_cart2(N12_coef,c.s.sample_rate);
  N13 = sample_state_cart2(N13_coef,c.s.sample_rate);
  N22 = sample_state_cart2(N22_coef,c.s.sample_rate);
  N23 = sample_state_cart2(N23_coef,c.s.sample_rate);
  N33 = sample_state_cart2(N33_coef,c.s.sample_rate);
  P11 = N11 - K11;
  P12 = N12 - K12;
  P13 = N13 - K13;
  P22 = N22 - K22;
  P23 = N23 - K23;
  P33 = N33 - K33;
end

function [u1, u2, u3, T11, T12, T13, T22, T23, T33, n_inv] ...
  = get_TT(c, state_coef,species);

  [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
    = get_primVars10(c, state_coef,species);
  ms = get_species_mass(c, species);
  n_inv = ms./rho;
  T11 = P11.*n_inv;
  T12 = P12.*n_inv;
  T13 = P13.*n_inv;
  T22 = P22.*n_inv;
  T23 = P23.*n_inv;
  T33 = P33.*n_inv;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sampling routines %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function data=sample_E1(c, state_coef); data=get_EN(c, state_coef, 1); end
function data=sample_E2(c, state_coef); data=get_EN(c, state_coef, 2); end
function data=sample_E3(c, state_coef); data=get_EN(c, state_coef, 3); end

function data = get_EN(c, state_coef, num)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  if(is_symmetric_pair_model(c))
    assert(num==3);
    coef = state_coef(:,:,:,stateIndices_s.E);
  else
    coef = state_coef(:,:,:,stateIndices_s.E(num));
  end
  data = sample_state_cart2(coef,c.s.sample_rate);
end

function [B1, B2, B3] = get_B(c, state_coef)
  B1 = sample_B1(c, state_coef);
  B2 = sample_B2(c, state_coef);
  B3 = sample_B3(c, state_coef);
end

function data=sample_B1(c, state_coef); data=get_BN(c, state_coef, 1); end
function data=sample_B2(c, state_coef); data=get_BN(c, state_coef, 2); end
function data=sample_B3(c, state_coef); data=get_BN(c, state_coef, 3); end

function data = get_BN(c, state_coef, num)
  if(num==3 && is_symmetric_pair_model(c))
    data = sample_zeros(c,1);
    return;
  end
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  coef = state_coef(:,:,:,stateIndices_s.B(num));
  data = sample_state_cart2(coef,c.s.sample_rate);
end

function data = sample_Eck(c, state_coef)
  % There is a lot of redundant computation here
  data = sample_Bxu(c, state_coef) ...
       + sample_pek(c, state_coef) ...
       + sample_dtJ(c, state_coef) ...
       - sample_E  (c, state_coef);
end

% Here we just assume that the inertial term makes up the difference
% in the momentum equation.
function data = sample_dtui(c, state_coef)
  % ideally should compare this with time differences
  % to confirm and test for numerical resistivity.
  % It would be good to demonstrate that numerical
  % resistivity is indeed anomalous.

  % There is redundant computation here:
  % we calculate the momentum and density
  pek = sample_pek(c, state_coef); % both here
  Bxu = sample_Bxu(c, state_coef); % and here.
  E = sample_state_var_from_state(c, 'E', state_coef);
  data = E - Bxu - pek;
end

function data = sample_JdotEprime(c, state_coef)
  nonideal = sample_nonideal(c,state_coef);
  J = sample_J(c,state_coef);
  assert(all(size(J)==size(nonideal)));
  data = zeros(size(J,1),size(J,2),1);
  for i=1:size(J,3)
    data = data + J(:,:,i).*nonideal(:,:,i);
  end
end

function data = sample_nonideal(c,state_coef)
  Bxu = sample_Bxu(c, state_coef);
  E = sample_state_var_from_state(c, 'E', state_coef);
  data = E - Bxu;
end

function data = sample_resistivity(c,state_coef)
  nonideal = sample_nonideal(c,state_coef);
  J = sample_J(c,state_coef);
  data = nonideal./J;
end

function data = sample_J3(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params = c.s.params;
  assert(isfield(stateIndices_s,'Mi'))
  if(isfield(stateIndices_s,'Me'))
    J_coef = state_coef(:,:,:,stateIndices_s.Mi(3))/params.ion_mass ...
            - state_coef(:,:,:,stateIndices_s.Me(3))/params.elc_mass;
  else
    J_coef = 2*state_coef(:,:,:,stateIndices_s.Mi(3))/params.ion_mass;
  end
  data = sample_state_cart2(J_coef,c.s.sample_rate);
end

function data = sample_J(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params = c.s.params;
  % need to implement J=curl(B)/mu_0 for mhd model
  assert(isfield(stateIndices_s,'Mi'));
  if(is_symmetric_pair_model(c))
    J_coef = 2*state_coef(:,:,:,stateIndices_s.Mi(3))/params.ion_mass;
  else
    assert(isfield(stateIndices_s,'Me'))
    J_coef = state_coef(:,:,:,stateIndices_s.Mi)/params.ion_mass ...
            - state_coef(:,:,:,stateIndices_s.Me)/params.elc_mass;
  end
  data = sample_state_cart2(J_coef,c.s.sample_rate);
end

function [Ji, Je, J] = sample_Jie(c, state_coef)
  params = c.s.params;
  ion_mass = params.ion_mass;
  elc_mass = params.elc_mass;
  [ui, rho_i, Mi] = sample_ui(c, state_coef);
  [ue, rho_e, Me] = sample_ue(c, state_coef);
  M = Mi + Me;
  rho = rho_i + rho_e;
  u = tupleOverScalar(M, rho);
  J = Mi/ion_mass - Me/elc_mass;
  %
  % there are three natural reference frames
  % in which we might wish to know the species currents:
  % (1) the frame of the grid,
  % (2) the frame of the fluid, and
  % (3) the frame of the flux-transporting flow.
  %
  % compute the reference frame in which an ideal
  % magnetic field would be convected according to
  % Ohm's law with the Hall term
  %
  u0 = u + tupleOverScalar(J,rho./(elc_mass-ion_mass));
  wi = ui - u0;
  % charge density
  sigma_i = rho_i/ion_mass;
  Ji = tupleTimesScalar(wi, sigma_i);
  %
  we = ue - u0;
  sigma_e = rho_e/(-elc_mass);
  Je = tupleTimesScalar(we, sigma_e);
end

function Ji = sample_Ji(c, state_coef); [Ji,Je]=sample_Jie(c, state_coef); end
function Je = sample_Je(c, state_coef); [Ji,Je]=sample_Jie(c, state_coef); end

function [Mi, Me] = sample_Mie(c, state_coef);
  Mi = sample_state_var_from_state(c, 'Mi', state_coef);
  Me = sample_state_var_from_state(c, 'Me', state_coef);
end

function [rho_i, rho_e] = sample_rho_ie(c, state_coef);
  rho_i = sample_state_var_from_state(c, 'rho_i', state_coef);
  rho_e = sample_state_var_from_state(c, 'rho_e', state_coef);
end

function xi = sample_xi(c, state_coef)
  rho_i = sample_state_var_from_state(c, 'rho_i', state_coef);
  x_i = sample_state_var_from_state(c, 'x_i', state_coef);
  xi = x_i./rho_i;
end

function yi = sample_yi(c, state_coef)
  rho_i = sample_state_var_from_state(c, 'rho_i', state_coef);
  y_i = sample_state_var_from_state(c, 'y_i', state_coef);
  yi = y_i./rho_i;
end

function vals = sample_Mi1(c, state_coef); vals = sample_MsN(c, state_coef, 'i', 1); end
function vals = sample_Mi2(c, state_coef); vals = sample_MsN(c, state_coef, 'i', 2); end
function vals = sample_Mi3(c, state_coef); vals = sample_MsN(c, state_coef, 'i', 3); end
function vals = sample_Me1(c, state_coef); vals = sample_MsN(c, state_coef, 'e', 1); end
function vals = sample_Me2(c, state_coef); vals = sample_MsN(c, state_coef, 'e', 2); end
function vals = sample_Me3(c, state_coef); vals = sample_MsN(c, state_coef, 'e', 3); end

function MsN_vals = sample_MsN(c, state_coef, species, component)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  space_order = c.s.params.space_order;
  Ms = ['M' species];
  MsN_coef = state_coef(:,:,:,stateIndices_s.(Ms)(component));
  MsN_vals = sample_state_cart2(MsN_coef,c.s.sample_rate);
end

function vals = sample_ui1(c, coef); vals = sample_usN(c, coef, 'i', 1); end
function vals = sample_ui2(c, coef); vals = sample_usN(c, coef, 'i', 2); end
function vals = sample_ui3(c, coef); vals = sample_usN(c, coef, 'i', 3); end
function vals = sample_ue1(c, coef); vals = sample_usN(c, coef, 'e', 1); end
function vals = sample_ue2(c, coef); vals = sample_usN(c, coef, 'e', 2); end
function vals = sample_ue3(c, coef); vals = sample_usN(c, coef, 'e', 3); end

function vals = sample_usN(c, state_coef, species, component, sample_rate);
  if(~exist('sample_rate','var')); sample_rate = c.s.sample_rate; end;
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  space_order = c.s.params.space_order;
  Ms = ['M' species];
  MsN_coef = state_coef(:,:,:,stateIndices_s.(Ms)(component));
  MsN_vals = sample_state_cart2(MsN_coef,sample_rate);
  rho_s_vals = sample_rho_s(c, state_coef, species, sample_rate);
  vals = MsN_vals./rho_s_vals;
end

function [ui, rho_i, Mi] = sample_ui(c, state_coef);
  [ui, rho_i, Mi] = sample_u_s(c, state_coef, 'i');
end
function [ue, rho_e, Me] = sample_ue(c, state_coef);
  [ue, rho_e, Me] = sample_u_s(c, state_coef, 'e');
end

function [u_s, rho_s, M_s] = sample_u_s(c, state_coef, species);
  str_Ms = ['M' species];
  str_rho_s = ['rho_' species];
  M_s  = sample_state_var_from_state(c, str_Ms, state_coef);
  rho_s  = sample_state_var_from_state(c, str_rho_s, state_coef);
  numdims = size(M_s,3);
  u_s = zeros(size(M_s));
  for i=1:numdims; u_s(:,:,i) = M_s(:,:,i)./rho_s; end
end

function val = sample_Ms(c, state_coef, species)
  val = sample_state_var_from_state(c, ['M' species], state_coef);
end

function val = sample_Mi(c, state_coef)
  val = sample_Ms(c, state_coef, 'i');
end

function val = sample_Me(c, state_coef)
  if(is_symmetric_pair_model(c))
    val = sample_Ms(c, state_coef, 'i');
    val(:,:,3) = - val(:,:,3);
  else
    val = sample_Ms(c, state_coef, 'e');
  end
end

function M_val = sample_M(c, state_coef)
  indices_s = c.s.stateVarMap.stateIndices_s;
  if(is_symmetric_pair_model(c))
    indices = indices_s.Mi(1:2); % third component is zero.
    M_state = 2.*state_coef(:,:,:,indices);
  else
    M_state = state_coef(:,:,:,indices_s.Mi) ...
            + state_coef(:,:,:,indices_s.Me);
  end
  M_val = sample_state_cart2(M_state,c.s.sample_rate);
end

function rho_val = sample_rho(c, state_coef)
  indices_s = c.s.stateVarMap.stateIndices_s;
  if(is_symmetric_pair_model(c))
    rho_state = 2.*state_coef(:,:,:,indices_s.rho_i);
  else
    rho_state = state_coef(:,:,:,indices_s.rho_i) ...
              + state_coef(:,:,:,indices_s.rho_e);
  end
  rho_val = sample_state_cart2(rho_state,c.s.sample_rate);
end

function [u, rho, M] = sample_u(c, state_coef);
  M = sample_M(c, state_coef);
  rho = sample_rho(c, state_coef);
  numdims = size(M,3);
  u = zeros(size(M));
  for j=1:numdims
    u(:,:,j) = M(:,:,j)./rho;
  end
end

function vals = sample_ni(c, state_coef); vals = get_ns(c, state_coef, 'i'); end
function vals = sample_ne(c, state_coef); vals = get_ns(c, state_coef, 'e'); end

function vals = get_ns(c, state_coef, species);
  rho_vals = sample_rho_s(c, state_coef, species);
  m_s = get_species_mass(c, species);
  vals = rho_vals./m_s;
end

function rho_vals = sample_rho_s(c, state_coef, species, sample_rate);
  if(~exist('sample_rate','var')); sample_rate = c.s.sample_rate; end;
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  space_order = c.s.params.space_order;
  rho = ['rho_' species];
  rho_coef = state_coef(:,:,:,stateIndices_s.(rho));
  rho_vals  = sample_state_cart2(rho_coef,sample_rate);
end

function u_val = sample_u_val(M_state,rho_state, how_sample);
  numdims = size(M_state,4);
  M_val = sample_state_cart2(M_state,how_sample);
  rho_val = sample_state_cart2(rho_state,how_sample);
  u_val = zeros(size(M_val));
  for j=1:numdims
    u_val(:,:,j) = M_val(:,:,j)./rho_val;
  end
end

% sample the Hall term
function data = sample_hall(c, state_coef);
  params = c.s.params;
  JxB = sample_JxB(c, state_coef);
  rho = sample_rho(c, state_coef);
  data = tupleOverScalar((params.ion_mass - params.elc_mass).*JxB,rho);
end

function JxB = sample_JxB(c, state_coef);
  J = sample_J(c, state_coef);
  B = sample_state_var_from_state(c, 'B', state_coef);
  JxB = val_crossProduct(J,B);
end

function [Bxus, rho_s] = sample_Bxus(c, state_coef, species);
  [u_s, rho_s, Ms] = sample_u_s(c, state_coef, species);
  B = sample_state_var_from_state(c, 'B', state_coef);
  Bxus = val_crossProduct(B,ue);
end

function [Bxui, rho_i] = sample_Bxui(c, state_coef);
  Bxui = sample_Bxus(c, state_coef, 'i');
end

function [Bxue, rho_e] = sample_Bxue(c, state_coef);
  Bxue = sample_Bxus(c, state_coef, 'e');
end

function out = sample_Bxu3(c, state_coef);
  u = sample_u(c, state_coef);
  B = sample_state_var_from_state(c, 'B', state_coef);
  out = val_crossProduct3(B,u);
end

function [Bxus3, rho_s] = sample_Bxus3(c, state_coef, species);
  [u_s, rho_s, Ms] = sample_u_s(c, state_coef, species);
  B = sample_state_var_from_state(c, 'B', state_coef);
  Bxus3 = val_crossProduct3(B,u_s);
end

function [Bxui3, rho_i] = sample_Bxui3(c, state_coef);
  [Bxui3, rho_i] = sample_Bxus3(c, state_coef, 'i');
end

function Bxue3 = sample_Bxue3(c, state_coef);
  [Bxue3, rho_e] = sample_Bxus3(c, state_coef, 'e');
end

% this currently is only implemented for symmetric pair plasma
%
function data = sample_Bxu(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  B_coef = state_coef(:,:,:,stateIndices_s.B);
  B_val = sample_state_cart2(B_coef,c.s.sample_rate);
  if(is_symmetric_pair_model(c))
    assert(numel(stateIndices_s.B)==2)
    % returns a scalar
    % magnetic field lies in the plane
    assert(numel(stateIndices_s.B)==2);
    rho_i_coef = state_coef(:,:,:,stateIndices_s.rho_i);
    Mi_coef = state_coef(:,:,:,stateIndices_s.Mi([1 2]));
    u_val = sample_u_val(Mi_coef, rho_i_coef, c.s.sample_rate);
    data = val_crossProduct3(B_val,u_val);
  else
    assert(numel(stateIndices_s.B)==3)
    rho_i_coef = state_coef(:,:,:,stateIndices_s.rho_i);
    rho_e_coef = state_coef(:,:,:,stateIndices_s.rho_e);
    rho_coef = rho_i_coef + rho_e_coef;
    Mi_coef = state_coef(:,:,:,stateIndices_s.Mi);
    Me_coef = state_coef(:,:,:,stateIndices_s.Me);
    M_coef = Mi_coef + Me_coef;
    u_val = sample_u_val(M_coef, rho_coef, c.s.sample_rate);
    data = val_crossProduct(B_val,u_val);
  end
end

% I need to make a version for symmetric pair plasma
%
function curlB = sample_curlB(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params = c.s.params;
  space_order = params.space_order;
  cBx = state_coef(:,:,:,stateIndices_s.B(1));
  cBy = state_coef(:,:,:,stateIndices_s.B(2));
  cBz = state_coef(:,:,:,stateIndices_s.B(3));
  % do I have the flip symmetries right here?
  cDyBz = dy_state(c,cBz,space_order,params.plot_dy,[1,1]);
  cDyBx = dy_state(c,cBx,space_order,params.plot_dy,[-1,1]); % -1,1 works
  cDxBz = dx_state(c,cBz,space_order,params.plot_dy,[1,1]);
  cDxBy = dx_state(c,cBy,space_order,params.plot_dy,[-1,-1]);
  DyBz = sample_state_cart2( cDyBz, c.s.sample_rate);
  curlBz = sample_state_cart2( cDxBy - cDyBx, c.s.sample_rate);
  curlB = zeros(size(curlBz,1),size(curlBz,2),3);
  curlB(:,:,1) = DyBz;
  curlB(:,:,2) = sample_state_cart2(-cDxBz, c.s.sample_rate);
  curlB(:,:,3) = curlBz;
end

function out = sample_Pi13x(c, coef); out = sample_Ps13x(c, coef, 'i'); end
function out = sample_Pe13x(c, coef); out = sample_Ps13x(c, coef, 'e'); end

function Ps13x_vals = sample_Ps13x(c, state_coef, species)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params=c.s.params;
  space_order = params.space_order;
  %
  % get kinetic energy state coefficients
  %
  %[rho_s_, Ms_] = get_gasVarNames4(species);
  [rho_s_, Ms_, Ns_] = get_gasVarNames10(species);
  Ms1_coef = state_coef(:,:,:,stateIndices_s.(Ms_)(1));
  Ms3_coef = state_coef(:,:,:,stateIndices_s.(Ms_)(3));
  Ms1_qval = sample_state_cart2(Ms1_coef,0);
  Ms3_qval = sample_state_cart2(Ms3_coef,0);
  rho_s_coef = state_coef(:,:,:,stateIndices_s.(rho_s_));
  rho_s_qval  = sample_state_cart2(rho_s_coef,0);
  us3_qval = Ms3_qval./rho_s_qval;
  Ns13_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(3));
  %Ns12_3_coef = state_coef(:,:,:,stateIndices_s.Ns([3 5]));
  Ks13_qval = Ms1_qval.*us3_qval;
  Ks13_coef = project_onto_basis_cart2(Ks13_qval,space_order);
  %
  % get ion pressure state coefficients
  Ps13_coef = Ns13_coef - Ks13_coef;
  Ps13x_coef = dx_state(c,Ps13_coef,space_order,params.plot_dx,[-1,-1]);
  Ps13x_vals = sample_state_cart2(Ps13x_coef,c.s.sample_rate);
end

function out = sample_Pi23y(c, coef); out = sample_Ps23y(c, coef, 'i'); end
function out = sample_Pe23y(c, coef); out = sample_Ps23y(c, coef, 'e'); end

function Ps23y_vals = sample_Ps23y(c, state_coef, species)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params=c.s.params;
  space_order = params.space_order;
  %
  % get kinetic energy state coefficients
  %
  [rho_s_, Ms_, Ns_] = get_gasVarNames10(species);
  Ms2_coef = state_coef(:,:,:,stateIndices_s.(Ms_)(2));
  Ms3_coef = state_coef(:,:,:,stateIndices_s.(Ms_)(3));
  Ms2_qval = sample_state_cart2(Ms2_coef,0);
  Ms3_qval = sample_state_cart2(Ms3_coef,0);
  rho_s_coef = state_coef(:,:,:,stateIndices_s.(rho_s_));
  rho_s_qval  = sample_state_cart2(rho_s_coef,0);
  us3_qval = Ms3_qval./rho_s_qval;
  Ns23_coef = state_coef(:,:,:,stateIndices_s.(Ns_)(5));
  Ks23_qval = Ms2_qval.*us3_qval;
  Ks23_coef = project_onto_basis_cart2(Ks23_qval,space_order);
  %
  % get ion pressure state coefficients
  Ps23_coef = Ns23_coef - Ks23_coef;
  Ps23y_coef = dy_state(c,Ps23_coef,space_order,params.plot_dy,[-1,-1]);
  Ps23y_vals = sample_state_cart2(Ps23y_coef,c.s.sample_rate);
end

% compute the third component of the divergence of the
% ion pressure term divided by the ion charge density
%
function data_vals = sample_pek(c, state_coef, species)
  divPs3_vals = sample_divPs3(c, state_coef, species);
  rho_s_ = get_gasVarNames4(species);
  rho_s_coef = state_coef(:,:,:,stateIndices_s.(rho_s_));
  rho_s_vals = sample_state_cart2(rho_s_coef);
  ns = rho_s_vals./m_s;
  data_vals = tupleOverScalar(divPi3_vals,ns_vals);
end

function out = sample_divPi3(c, state_coef)
  out = sample_divPs3(c, state_coef, 'i');
end
function out = sample_divPe3(c, state_coef)
  out = sample_divPs3(c, state_coef, 'e');
end

function divB_vals = sample_divB(c, state_coef)
  B1x_vals = sample_B1x(c, state_coef);
  B2y_vals = sample_B2y(c, state_coef);
  divB_vals = B1x_vals + B2y_vals;
end

function B1x_vals = sample_B1x(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params=c.s.params;
  space_order = params.space_order;

  B1_coef = state_coef(:,:,:,stateIndices_s.B(1));
  B1x_coef = dx_state(c,B1_coef,space_order,params.plot_dx,[1,1]);
  B1x_vals = sample_state_cart2(B1x_coef,c.s.sample_rate);
end

function B2y_vals = sample_B2y(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params=c.s.params;
  space_order = params.space_order;

  B2_coef = state_coef(:,:,:,stateIndices_s.B(2));
  B2y_coef = dy_state(c,B2_coef,space_order,params.plot_dy,[1,-1]);
  B2y_vals = sample_state_cart2(B2y_coef,c.s.sample_rate);
end

function divEc = sample_divEc(c, state_coef)
  divE = sample_divE(c, state_coef);
  ni = sample_ni(c, state_coef);
  ne = sample_ne(c, state_coef);
  sigma = ni-ne;
  divEc = divE - sigma*c.s.params.one_over_epsilon;
end

function divE_vals = sample_divE(c, state_coef)
  E1x_vals = sample_E1x(c, state_coef);
  E2y_vals = sample_E2y(c, state_coef);
  divE_vals = E1x_vals + E2y_vals;
end

function E1x_vals = sample_E1x(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params=c.s.params;
  space_order = params.space_order;

  E1_coef = state_coef(:,:,:,stateIndices_s.E(1));
  E1x_coef = dx_state(c,E1_coef,space_order,params.plot_dx,[-1,-1]);
  E1x_vals = sample_state_cart2(E1x_coef,c.s.sample_rate);
end

function E2y_vals = sample_E2y(c, state_coef)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params=c.s.params;
  space_order = params.space_order;

  E2_coef = state_coef(:,:,:,stateIndices_s.E(2));
  E2y_coef = dy_state(c,E2_coef,space_order,params.plot_dy,[-1,-1]);
  E2y_vals = sample_state_cart2(E2y_coef,c.s.sample_rate);
end

function divPs3_vals = sample_divPs3(c, state_coef, species)
  alternate_implementation=0;
  % this involves a redundant calculation of Ms1_qval
  if(alternate_implementation)
    divPs3_vals = sample_Ps13x(c, state_coef, species) ...
                + sample_Ps23y(c, state_coef, species);
    return;
  end

  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  params=c.s.params;
  space_order = params.space_order;
  [rho_s_, Ms_, Ns_] = get_gasVarNames10(species);
  %
  % get kinetic energy state coefficients
  %
  Ms_coef = state_coef(:,:,:,stateIndices_s.(Ms_));
  Ms_qval = sample_state_cart2(Ms_coef,0);
  rho_s_coef = state_coef(:,:,:,stateIndices_s.(rho_s_));
  rho_s_qval = sample_state_cart2(rho_s_coef,0);

  us3_qval = Ms_qval(:,:,3)./rho_s_qval;
  sz=size(Ms_qval); sz(3)=2;
  Ks12_3_qval = zeros(sz);
  Ks12_3_qval(:,:,1) = Ms_qval(:,:,1).*us3_qval;
  Ks12_3_qval(:,:,2) = Ms_qval(:,:,2).*us3_qval;
  Ks12_3_coef = project_onto_basis_cart2(Ks12_3_qval,space_order);
  %
  % get ion pressure state coefficients
  Ns12_3_coef = state_coef(:,:,:,stateIndices_s.(Ns_)([3 5]));
  Ps12_3_coef = Ns12_3_coef - Ks12_3_coef;
  % take the divergence of the pressure
  divPs3_coef ...
    = dx_state(c,Ps12_3_coef(:,:,:,1),space_order,params.plot_dx,[-1,-1]) ... % Ps13x
    + dy_state(c,Ps12_3_coef(:,:,:,2),space_order,params.plot_dy,[-1,-1]);    % Ps23y
  divPs3_vals=sample_state_cart2(divPs3_coef,c.s.sample_rate);
  divPs3_vals(:,end)=0; % suppress edge effects
  m_s = get_species_mass(c, species);
  %rho_s_vals = sample_state_cart2(rho_s_coef);
  %ns = rho_s_vals./m_s;
end

function out = sample_divui(c, coef); out = sample_divus(c, coef, 'i'); end
function out = sample_divue(c, coef); out = sample_divus(c, coef, 'e'); end

function divus = sample_divus(c, state_coef, species)
  us1x = sample_us1x(c, state_coef, species);
  us2y = sample_us2y(c, state_coef, species);
  divus = us1x + us2y;
end

function out = sample_ui1x(c, coef); out = sample_us1x(c, coef, 'i'); end
function out = sample_ui2y(c, coef); out = sample_us2y(c, coef, 'i'); end
function out = sample_ue1x(c, coef); out = sample_us1x(c, coef, 'e'); end
function out = sample_ue2y(c, coef); out = sample_us2y(c, coef, 'e'); end

function us2y = sample_us2y(c, state_coef, species)
  % get u3 at quadrature points
  u2s_qval = sample_usN(c, state_coef, species, 2, 0);
  % project onto legendre basis
  params = c.s.params;
  space_order = params.space_order;
  u2s_coef = project_onto_basis_cart2(u2s_qval,space_order);
  % differentiate
  u2s_y_coef = dy_state(c,u2s_coef,space_order,params.plot_dy,[-1,-1]);
  % sample at plotting points
  us2y = sample_state_cart2(u2s_y_coef, c.s.sample_rate);
end

function us1x = sample_us1x(c, state_coef, species)
  % get u3 at quadrature points
  u1s_qval = sample_usN(c, state_coef, species, 1, 0);
  % project onto legendre basis
  params = c.s.params;
  space_order = params.space_order;
  u1s_coef = project_onto_basis_cart2(u1s_qval,space_order);
  % differentiate
  u1s_x_coef = dx_state(c,u1s_coef,space_order,params.plot_dx,[-1,1]);
  % sample at plotting points
  us1x = sample_state_cart2(u1s_x_coef, c.s.sample_rate);
end

% u_s.Grad(u3_s)
function udotDu3 = sample_udotDu3_s(c, state_coef, species)
  % get gradient of u3
  %
  % get u3 at quadrature points
  u3s_qval = sample_usN(c, state_coef, species, 3, 0);
  % project onto legendre basis
  params = c.s.params;
  space_order = params.space_order;
  u3s_coef = project_onto_basis_cart2(u3s_qval,space_order);
  % differentiate
  u3s_x_coef = dx_state(c,u3s_coef,space_order,params.plot_dx,[1,1]);
  u3s_y_coef = dy_state(c,u3s_coef,space_order,params.plot_dy,[1,1]);
  % sample at plotting points
  u3s_x_vals = sample_state_cart2(u3s_x_coef, c.s.sample_rate);
  u3s_y_vals = sample_state_cart2(u3s_y_coef, c.s.sample_rate);

  % sample u at plotting points
  u1s_qval = sample_usN(c, state_coef, species, 1);
  u2s_qval = sample_usN(c, state_coef, species, 2);
  % compute u.Grad u at plotting points
  udotDu3 = u1s_qval.*u3s_x_vals + ...
            u2s_qval.*u3s_y_vals;
end

function out = sample_udotDu3_i(c, state_coef)
  out=sample_udotDu3_s(c, state_coef, 'i');
end
function out = sample_udotDu3_e(c, state_coef)
  out=sample_udotDu3_s(c, state_coef, 'e');
end

function p = sample_p1i(c, state_coef); p = get_pNs(c, state_coef, 'i', 1); end
function p = sample_p2i(c, state_coef); p = get_pNs(c, state_coef, 'i', 2); end
function p = sample_p3i(c, state_coef); p = get_pNs(c, state_coef, 'i', 3); end
function p = sample_p1e(c, state_coef); p = get_pNs(c, state_coef, 'e', 1); end
function p = sample_p2e(c, state_coef); p = get_pNs(c, state_coef, 'e', 2); end
function p = sample_p3e(c, state_coef); p = get_pNs(c, state_coef, 'e', 3); end

function T = sample_T1i(c, state_coef); T = get_TNs(c, state_coef, 'i', 1); end
function T = sample_T2i(c, state_coef); T = get_TNs(c, state_coef, 'i', 2); end
function T = sample_T3i(c, state_coef); T = get_TNs(c, state_coef, 'i', 3); end
function T = sample_T1e(c, state_coef); T = get_TNs(c, state_coef, 'e', 1); end
function T = sample_T2e(c, state_coef); T = get_TNs(c, state_coef, 'e', 2); end
function T = sample_T3e(c, state_coef); T = get_TNs(c, state_coef, 'e', 3); end

function TN = get_TNs(c, state_coef, species, num);
  [p, rho] = get_pNs(c, state_coef, species, num);
  m_s = get_species_mass(c, species);
  TN = m_s*p./rho;
end

function [p, rho] = get_pNs(c, state_coef, species, num);
  [p1, p2, p3, rho] = sample_eigs_s(c, state_coef, species);
  switch(num)
    case 1; p=p1;
    case 2; p=p2;
    case 3; p=p3;
    otherwise; error('invalid value: num = ', num2str(num));
  end
end

function bPb = sample_pbi(c, state); bPb = get_pbs(c, state, 'i'); end
function bPb = sample_pbe(c, state); bPb = get_pbs(c, state, 'e'); end
function bPb = get_pbs(c, state_coef, species)
  [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
     = get_primVars10(c, state_coef, species);
  [B1, B2, B3] = get_B(c, state_coef);
  %
  B12 = B1.*B2;
  B13 = B1.*B3;
  B23 = B2.*B3;
  B11 = B1.^2;
  B22 = B2.^2;
  B33 = B3.^2;
  B_mag2 = B11 + B22 + B33;
  one_over_Bmag2 = 1./B_mag2;
  %
  BPB = B11.*P11 + B22.*P22 + B33.*P33 + 2*(B12.*P12 + B13.*P13 + B23.*P23);
  bPb = BPB.*one_over_Bmag2; % b.P.b
end

% could use this for pressure or temperature
%
% This is not enough as a test of gyrotropy
%
function [para, perp1, perp2] ...
  = get_gyrotropic_components(B1, B2, B3, P11, P12, P13, P22, P23, P33)

  % project onto B
  %
  %B12 = B1.*B2;
  %B13 = B1.*B3;
  %B23 = B2.*B3;
  B11 = B1.^2;
  B22 = B2.^2;
  B33 = B3.^2;
  B_mag2 = B11 + B22 + B33;
  one_over_Bmag2 = 1./B_mag2;
  % These aliases incur no cost since matlab uses copy-on-write
  % and make formulas easier to enter.
  P21 = P12; P31 = P13; P32 = P23;
  PB_1 = P11.*B1 + P12.*B2 + P13.*B3; % (P.B)_1
  PB_2 = P21.*B1 + P22.*B2 + P23.*B3; % (P.B)_2
  PB_3 = P31.*B1 + P32.*B2 + P33.*B3; % (P.B)_3
  BPPB = PB_1.*PB_1 + PB_2.*PB_2 + PB_3.*PB_3; % B.P.P.B
  BPB = B1.*PB_1 + B2.*PB_2 + B3.*PB_3; % B.P.B
  %BPB = B11.*P11 + B22.*P22 + B33.*P33 + 2*(B12.*P12 + B13.*P13 + B23.*P23);
  bPPb = BPPB.*one_over_Bmag2; % b.P.P.b
  bPb = BPB.*one_over_Bmag2; % b.P.b
  trP = P11+P22+P33; % trace(P)
  trP_2 = trP*.5;
  PP_2 = 0.5*(P11.*P11+P22.*P22+P33.*P33) + P12.*P12+P13.*P13+P23.*P23; % P:P
  %
  % find extrema of n.P.n such that n.B = 0.
  % Theory:
  %   Let k = b = Bhat.
  %   Let i,j,k be an orthonormal system.
  %   Let P11 = i.P.i, P12 = i.P.j, P22 = j.P.j
  %   Need the eigenvalues of the 2x2 matrix A:=P_{ij}.
  %   Need characteristic polynomial of A.
  %   Need trace and determinant of A.
  %   (want tensorial formulas without the somewhat arbitrary vectors i and j)
  %   tr(A) = (ii+jj):P = (I-kk):P = tr(P)-b.P.b
  trA = trP - bPb;
  %   For the determinant,
  %   recall the tensorial formula for a 2 by 2 determinant,
  %   2*det(A) = eps_{ij}*eps_{lm}*P_{il}*P_{jm}.
  %   To make this tensorial note that
  %   eps_{ij} = eps_{ijn} k_n = -k cross I.
  %   So we have
  %   2*det(A) = eps_{ijk}*k_k*eps_{lmn}*k_n*P_{il}*P_{jm}.
  %   (Note that in this formula it makes no difference if we allow i and j
  %   to range over the values (1,2,3) instead of merely (1,2).
  %   So we in fact have a tensorial expression for the determinant,
  %   2*det(A) = eps_{ijk}*b_k*eps_{lmn}*b_n*P_{il}*P_{jm}, i.e.,
  %   2*det(A) = eps_{ijk}*eps_{lmn}*P_{il}*P_{jm}*(bb)_{kn}
  %   which I remark looks very much like the formula
  %   3!*det(P) = eps_{ijk}*eps_{lmn}*P_{il}*P_{jm}*P_{kn}
  %   for the determinant of a 3x3 matrix.
  %   Recall that eps_{ijk}*eps_{lmn}
  %     = d_il*d_jm*d_kn - d_il*d_jn*d_km
  %     + d_im*d_jn*d_kl - d_im*d_jl*d_kn
  %     + d_in*d_jl*d_km - d_in*d_jm*d_kl.
  %   where we have abbreviated d_ij = delta_{ij}, etc. So
  %   So making the definition b_ij = b_i*b_j,
  %   2*det(A)   
  %    =  P_ii*P_jj*b_kk - P_ii*P_jk*b_kj
  %     + P_ik*P_ji*b_kj - P_ik*P_jj*b_ki
  %     + P_ij*P_jk*b_ki - P_ij*P_ji*b_kk
  %    =  (tr P)^2 - (tr P) b.P.b
  %      + b.P.P.b - (tr P) b.P.b
  %      + b.P.P.b - P:P
  %    =  (tr P)^2 - 2*(tr P)*P_bb
  %      + 2*(b.P.P.b) - P:P
  % So
  %   det(A) =  (tr P)*((tr P)/2 - P_bb)
  %             + (b.P.P.b) - P:P/2
  detA = trP.*(trP_2 - bPb) + bPPb - PP_2;
  % characteristic polynomial is
  %   lambda^2 - 2*(tr(A)/2)+det(A)
  % roots are
  %   lambda = tr(A)/2 +/- sqrt((tr(A)/2)^2-det(A))
  trA_2 = trA.*.5;
  discriminant = trA_2.*trA_2 - detA;
  assert(all(discriminant>0));
  sqrt_discriminant = sqrt(discriminant);
  para = bPb;
  perp1 = trA_2 - sqrt_discriminant;
  perp2 = trA_2 + sqrt_discriminant;
end

function pmax = sample_pmaxe(c, state_coef)
  pmax = sample_pmaxs(c, state_coef, 'e');
end

function pmax = sample_pmaxi(c, state_coef)
  pmax = sample_pmaxs(c, state_coef, 'i');
end

% it seems that eig3 is always the max
%
function pmax = sample_pmaxs(c, state_coef, species)
  % is there a shortcut to getting the maximal eigenvalue?
  [eig1, eig2, eig3, rho] = sample_eigs_s(c, state_coef, species);
  % should probably also verify that all eigenvalues are positive.
  pmax = max(eig1,eig2);
  pmax = max(pmax,eig3);
end

function pmin = sample_pmine(c, state_coef)
  pmin = sample_pmins(c, state_coef, 'e');
end

function pmin = sample_pmini(c, state_coef)
  pmin = sample_pmins(c, state_coef, 'i');
end

% it seems that eig1 is always the max
%
function pmin = sample_pmins(c, state_coef, species)
  [eig1, eig2, eig3] = sample_eigs_s(c, state_coef,species);
  pmin = min(eig1,eig2);
  pmin = min(pmin,eig3);
end

function pratio = sample_pratio_i(c, state_coef)
  pratio = sample_pratio_s(c, state_coef, 'i');
end

function pratio = sample_pratio_e(c, state_coef)
  pratio = sample_pratio_s(c, state_coef, 'e');
end

function pratio = sample_pratio_s(c, state_coef, species)
  [eig1, eig2, eig3, rho, u1, u2, u3] ...
    = sample_eigs_s(c, state_coef, species);
  pmin = min(eig1,eig2);
  pmin = min(pmin,eig3);
  pmax = max(eig1,eig2);
  pmax = max(pmax,eig3);
  pratio=pmax./pmin;
end

function maxspd = sample_maxspd_i(c, state_coef)
  maxspd = sample_maxspd_s(c, state_coef, 'i');
end

function maxspd = sample_maxspd_e(c, state_coef)
  maxspd = sample_maxspd_s(c, state_coef, 'e');
end

function maxspd = sample_maxspd_s(c, state_coef, species)
  % is there a shortcut to getting the maximal eigenvalue?
  [cfe, rho, u1, u2, u3] = sample_cfs(c, state_coef, species);
  umag = sqrt(u1.*u1 + u2.*u2+ u3.*u3);
  % in the ten-moment case this is not actually the maximum speed;
  % it is an upper bound on the speed.
  maxspd = umag + cfe;
end

function [cf] = sample_cfi(c, state_coef)
  cf = sample_cfs(c, state_coef, 'i');
end

function [cf] = sample_cfe(c, state_coef)
  cf = sample_cfs(c, state_coef, 'e');
end

% sample fast gas speed
% need to find largest eigenvalue of electron pressure tensor
function [cf, rho, u1, u2, u3] = sample_cfs(c, state_coef, species)
  % is there a shortcut to getting the maximal eigenvalue?
  species_model=get_species_model(c, species);
  switch species_model
    case 10
      [eig1, eig2, eig3, rho, u1, u2, u3] ...
        = sample_eigs_s(c, state_coef, species);
      maxeig = max(eig1,eig2);
      maxeig = max(maxeig,eig3);
      cs = sqrt(maxeig./rho);
      cf = cs*sqrt(3.);
    case 5
      % get pressure and density
      [rho, u1, u2, u3, p] = get_primVars05(c, state_coef, species);
      cf=sqrt(5./3.)*sqrt(p./rho);
    otherwise
      error(['unsupported value: ' species_model]);
  end
end

function tau_s = sample_isorate_i(c, state); tau_s = get_isorate_s(c, state, 'i'); end
function tau_s = sample_isorate_e(c, state); tau_s = get_isorate_s(c, state, 'e'); end
function [isorate, p] = get_isorate_s(c, state_coef, species)
  m_s = get_species_mass(c, species);
  base_iso_period = c.s.params.base_iso_period;
  switch c.s.params.iso_period_type
  case 'constant'
    tau_s = get_base_iso_period(c, species);
    % multiply by an array of ones?
    tau_s = (1./tau_s)*sample_ones(c,1);
    p = get_p_s(c, state_coef, species);
    return;
  case 'det'
    [detT, p, n] = sample_detTs(c, state_coef, species);
  case 'trace'
    [T,n,p] = sample_T_s(c, state_coef, species);
    detT = T.^3;
  otherwise
    error(['unsupported iso_period_type: ' species]);
  end
  isorate = (n./sqrt(detT)).*(1./(sqrt(m_s)*base_iso_period));
end

function tau = sample_tau_i(c, state); tau = get_iso_period_s(c, state, 'i'); end
function tau = sample_tau_e(c, state); tau = get_iso_period_s(c, state, 'e'); end

function tau_s = sample_detPtau_i(c, state); tau_s = get_detPtau_s(c,state,'i'); end
function tau_s = sample_detPtau_e(c, state); tau_s = get_detPtau_s(c,state,'e'); end

% This returns the Alec (detP) isotropization period
function tau = get_detPtau_s(c, state_coef, species);
  %m_s = get_species_mass(c, species);
  %base_iso_period = c.s.params.base_iso_period;
  %assert(base_iso_period > 0);
  %[detT, p, n] = sample_detTs(c, state_coef, species);
  %tau = (sqrt(detT)./n).*(base_iso_period*sqrt(m_s));
  [tau_s, p] = get_iso_period_s(c, state_coef, species, 'det');
end

function tau = sample_ptau_i(c, state); tau = get_ptau_s(c, state, 'i'); end
function tau = sample_ptau_e(c, state); tau = get_ptau_s(c, state, 'e'); end
function tau = get_ptau_s(c, state_coef, species);
  [tau, p] = get_iso_period_s(c, state_coef, species, 'trace');
end

% assumes the Braginskii closure
function [tau_s, p] = get_iso_period_s(c, state_coef, species, type)
  m_s = get_species_mass(c, species);
  base_iso_period = c.s.params.base_iso_period;
  switch c.s.params.iso_period_type
  case 'constant'
    tau_s = get_base_iso_period(c, species);
    % multiply by an array of ones?
    tau_s = tau_s*sample_ones(c,1);
    p = get_p_s(c, state_coef, species);
    return;
  case 'det'
    [detT, p, n] = sample_detTs(c, state_coef, species);
  case 'trace'
    [T,n,p] = sample_T_s(c, state_coef, species);
    detT = T.^3;
  otherwise
    error(['unsupported iso_period_type: ' species]);
  end
  tau_s = (sqrt(detT)./n).*(sqrt(m_s)*base_iso_period);
end

function [mu] = sample_mu_i(c, state)
  [mu] = get_viscosity_s(c, state, 'i');
end
function [mu] = sample_mu_e(c, state)
  [mu] = get_viscosity_s(c, state, 'e');
end
function [mu_s] = get_viscosity_s(c, state_coef, species)
  %base_iso_period = c.s.params.base_iso_period;
  %assert(base_iso_period > 0);
  %m_s = get_species_mass(c, species);
  %T = sample_T_s(c, state_coef, species);
  %mu_s = (base_iso_period.*sqrt(m_s))*T.^(5/2);
  [iso_period, p] = get_iso_period_s(c, state_coef, species);
  mu_s = iso_period.*p;
end

function K = sample_K_i(c, state_coef)
  K = get_thermal_conductivity_s(c, state_coef, 'i');
end
function K = sample_K_e(c, state_coef)
  K = get_thermal_conductivity_s(c, state_coef, 'e');
end
function K_s = get_thermal_conductivity_s(c, state_coef, species)
  m_s = get_species_mass(c, species);
  %
  %base_iso_period = c.s.params.base_iso_period;
  %assert(base_iso_period > 0);
  %T = sample_T_s(c, state_coef, species);
  %K_s = (base_iso_period/sqrt(m_s))*T.^(5/2);
  %
  viscosity = get_viscosity_s(c, state_coef, species);
  K_s = viscosity.*(1/m_s);
end

function qs = sample_qs(c, state_coef, species)
  % sample temperature at quadrature points
  c0 = c;
  c0.s.sample_rate = 0; % cause sampling to be done at quadrature points
  [qv_T] = sample_T_s(c0, state_coef, species);
  % project temperature onto polynomial basis
  T_coef = project_onto_basis_cart2(qv_T,space_order);
  % compute gradient of temperature (correct end-symmetries?)
  Tx_coef = dx_state(c,T_coef,space_order,params.plot_dx,[-1,-1]);
  Ty_coef = dy_state(c,T_coef,space_order,params.plot_dx,[-1,-1]);
  % sample gradient at plotting points
  Tx = sample_state_cart2(Tx_coef,c.s.sample_rate);
  Ty = sample_state_cart2(Ty_coef,c.s.sample_rate);
  % multiply by thermal_conductivity
  K_s = get_thermal_conductivity_s(c, state_coef, species);
  %sz=size(Tx);
  %sz(3)=2;
  %qs = zeros(sz);
  qs = sample_zeros(c,2);
  qs(:,:,1) = -K_s.*Tx;
  qs(:,:,2) = -K_s.*Ty;
end

function T = sample_Ti(c, state_coef); T = sample_T_s(c, state_coef, 'i'); end
function T = sample_Te(c, state_coef); T = sample_T_s(c, state_coef, 'e'); end

function [T,n,p] = sample_T_s(c, state_coef, species);
  [rho, u1, u2, u3, p] = get_primVars05(c, state_coef, species);
  m_s = get_species_mass(c, species);
  n = rho*(1./m_s);
  T = p./n;
end

function entropy = sample_entropy_i(c, state_coef)
  entropy = sample_entropy_s(c, state_coef, 'i');
end

function entropy = sample_entropy_e(c, state_coef)
  entropy = sample_entropy_s(c, state_coef, 'e');
end

function Entropy = sample_Entropy_i(c, state_coef)
  Entropy = sample_Entropy_s(c, state_coef, 'i');
end

function Entropy = sample_Entropy_e(c, state_coef)
  Entropy = sample_Entropy_s(c, state_coef, 'e');
end

% compute the entropy density per volume
function Entropy = sample_Entropy_s(c, state_coef, species)
  [entropy, n] = sample_entropy_s(c, state_coef, species);
  Entropy = n.*entropy;
  totentropy = sum(sum(Entropy))*(c.s.params.plot_dx*c.s.params.plot_dx);
  disp(['total entropy for species ' species ' is ' num2str(totentropy)]);
end

function [TA11,TA12,TA13,TA22,TA23,TA33] = get_adjugate(T11,T12,T13,T22,T23,T33)
  TA11 = T22.*T33 - T23.*T23;
  TA12 = T23.*T13 - T12.*T33;
  TA13 = T12.*T23 - T22.*T13;
  TA22 = T11.*T33 - T13.*T13;
  TA23 = T12.*T13 - T11.*T23;
  TA33 = T11.*T22 - T12.*T12;
end

% compute determinant of symmetric tensor
function detT = get_determinant(T11, T12, T13, T22, T23, T33);
  % find a row of the adjugate
  TA11 = T22.*T33 - T23.*T23;
  TA12 = T23.*T13 - T12.*T33;
  TA13 = T12.*T23 - T22.*T13;
  % TA22 = T11.*T33 - T13.*T13;
  % TA23 = T12.*T13 - T11.*T23;
  % TA33 = T11.*T22 - T12.*T12;
  %
  % find the determinant
  detT = TA11.*T11 + TA12.*T12 + TA13.*T13;
end

function [data,badind]=forcePositive(data, factor);
    badind = find(data<=0);
  if(numel(badind>0))
    goodind = find(data>0);
    mingood = min(data(goodind));
    data(badind) = mingood*factor;
  end
end

% compute the entropy density per particle number
function [entropy, n] = sample_entropy_s(c, state_coef, species)
  species_model=get_species_model(c, species);
  switch species_model
  case 10
    [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
       = get_primVars10(c, state_coef, species);
    % compute the determinant of the pressure tensor.
    detP = get_determinant(P11, P12, P13, P22, P23, P33);
    %detP = forcePositive(detP);
    [detP,badind]=forcePositive(detP, exp(-2));
    if(numel(badind)>0)
      disp(['pressure tensor is not everywhere positive definite ' ...
            'for species ' species]);
    end
    log_detP = log(detP);
  case 5
    [rho, u1, u2, u3, p] = get_primVars05(c, state_coef, species);
    [p,badind]=forcePositive(p, exp(-2*3));
    if(numel(badind)>0)
      disp(['pressure is not everywhere positive' 'for species ' species]);
    end
    log_detP = 3*log(p);
  otherwise
    error(['unsupported model: ' species_model]);
  end
  ms = get_species_mass(c, species);
  n = rho/ms;
  entropy = 0.5*(log_detP-5*log(n));
end

function p = sample_p(c, state_coef)
  p_i = get_p_s(c, state_coef, 'i');
  p_e = get_p_s(c, state_coef, 'e');
  p = p_i + p_e;
end

function p = sample_pi(c, state_coef); p = get_p_s(c, state_coef, 'i'); end
function p = sample_pe(c, state_coef); p = get_p_s(c, state_coef, 'e'); end

function p = get_p_s(c, state_coef, species);
  [rho, u1, u2, u3, p] = get_primVars05(c, state_coef, species);
end

function val = sample_Ni11(c, state_coef); val = get_Nsij(c, state_coef, 'i',1,1); end
function val = sample_Ni22(c, state_coef); val = get_Nsij(c, state_coef, 'i',2,2); end
function val = sample_Ni33(c, state_coef); val = get_Nsij(c, state_coef, 'i',3,3); end
function val = sample_Ni12(c, state_coef); val = get_Nsij(c, state_coef, 'i',1,2); end
function val = sample_Ni13(c, state_coef); val = get_Nsij(c, state_coef, 'i',1,3); end
function val = sample_Ni23(c, state_coef); val = get_Nsij(c, state_coef, 'i',2,3); end
function val = sample_Ne11(c, state_coef); val = get_Nsij(c, state_coef, 'e',1,1); end
function val = sample_Ne22(c, state_coef); val = get_Nsij(c, state_coef, 'e',2,2); end
function val = sample_Ne33(c, state_coef); val = get_Nsij(c, state_coef, 'e',3,3); end
function val = sample_Ne12(c, state_coef); val = get_Nsij(c, state_coef, 'e',1,2); end
function val = sample_Ne13(c, state_coef); val = get_Nsij(c, state_coef, 'e',1,3); end
function val = sample_Ne23(c, state_coef); val = get_Nsij(c, state_coef, 'e',2,3); end

function Nsij = get_Nsij(c, state_coef, species,i,j)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  idx = c.s.tensorIdxMap.linearIdx(i,j);
  Ns_ = ['N' species];
  stateIdx=stateIndices_s.(Ns_)(idx);
  Nsij_coef = state_coef(:,:,:,stateIdx);
  Nsij = sample_state_cart2(Nsij_coef,c.s.sample_rate);
end

function val = get_Psij(c, state_coef, species,i,j)
  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  [rho_s_, Ms_, Ns_] = get_gasVarNames10(species);
  rho_coef = state_coef(:,:,:,stateIndices_s.(rho_s_));
  rho = sample_state_cart2(rho_coef,c.s.sample_rate);
  Mi_coef  = state_coef(:,:,:,stateIndices_s.(Ms_)(i));
  Mi = sample_state_cart2(Mi_coef,c.s.sample_rate);
  if(i==j)
    % Mj_coef = Mi_coef;
    Mj = Mi;
  else
    Mj_coef  = state_coef(:,:,:,stateIndices_s.(Ms_)(j));
    Mj = sample_state_cart2(Mj_coef,c.s.sample_rate);
  end
  Kij = Mi.*Mj./rho;
  idx = c.s.tensorIdxMap.linearIdx(i,j);
  stateIdx=stateIndices_s.(Ns_)(idx);
  Nij_coef = state_coef(:,:,:,stateIdx);
  Nij = sample_state_cart2(Nij_coef,c.s.sample_rate);
  val = Nij-Kij;
end

function [tau] = get_base_iso_period(c, species);
  iso_period_type = c.s.params.iso_period_type;
  switch iso_period_type
  case 'constant'
    switch species
    case 'i'
      tau = c.s.params.ion_iso_period;
    case 'e'
      tau = c.s.params.elc_iso_period;
    otherwise
      error(['unsupported species: ' species]);
    end
  case {'det','trace'}
    tau = c.s.params.base_iso_period;
  otherwise
    error(['unsupported iso_period_type: ' species]);
  end
end

function sit = sample_sit(c, state); sit = sample_sst(c, state, 'i'); end
function set = sample_set(c, state); set = sample_sst(c, state, 'e'); end
function sst = sample_sst(c, state_coef, species)
  % Let s = entropy per mass.
  %   2*rho*d_t s = TTinv:R, where R = (p*I-P)/tau.
  % But rho=n*m, p = n*T, and P = n*TT, so this says that
  %   2*m*d_t s = TTinv:(T*I-TT)/tau, i.e.
  %       d_t s = TTinv:(T*I-TT)/(tau*2*m).
  %             = (tr(TTinv)*tr(TT)-9)/(6*tau*m).
  [u1, u2, u3, T11, T12, T13, T22, T23, T33, n_inv] ...
    = get_TT(c,state_coef,species);
  [TA11,TA12,TA13,TA22,TA23,TA33] = get_adjugate(T11,T12,T13,T22,T23,T33);
  detT = get_determinant(T11, T12, T13, T22, T23, T33);
  trTTinv = (TA11+TA22+TA33)./detT;
  trTT = T11+T22+T33;
  m_s = get_species_mass(c, species);
  tau_s = get_tau_s(c, species, n_inv, detT, trTT);
  sst = (trTTinv.*trTT-9.)./(tau_s.*(6*m_s));
end

% a version of get_iso_period_s with different arguments
function tau_s = get_tau_s(c, species, n_inv, detT, trTT)
  m_s = get_species_mass(c, species);
  base_iso_period = c.s.params.base_iso_period;
  switch c.s.params.iso_period_type
  case 'constant'
    tau_s = get_base_iso_period(c, species);
    % multiply by an array of ones?
    tau_s = tau_s*sample_ones(c,1);
    return;
  case 'det'
    disp(['using detT with tau_0 = ' num2str(base_iso_period)]);
    % do nothing
  case 'trace'
    disp(['using trT with tau_0 = ' num2str(base_iso_period)]);
    T=trTT*(1./3.);
    detT = T.^3;
  otherwise
    error(['unsupported iso_period_type: ' species]);
  end
  tau_s = (sqrt(detT).*n_inv).*(sqrt(m_s)*base_iso_period);
end


function Sit = sample_Sit(c, state); Sit = sample_Sst(c, state, 'i'); end
function Set = sample_Set(c, state); Set = sample_Sst(c, state, 'e'); end
function Sst = sample_Sst(c, state_coef, species)
  if(0)
  sst = sample_sst(c, state_coef, species);
  rho_s = sample_rho_s(c, state_coef, species);
  Sst = rho_s.*sst;
  return;
  end
  % rate of entropy production is twice
  %
  % TTinv:(pI-P)/tau = (p*tr(TTinv)-3n)/tau
  % = (tr(TT)*tr(TTinv)-9)*n/(3*tau)
  % = (tr(P)*tr(Pinv)-9)*n/(3*tau).
  % If tau = tau_0*sqrt(m_s*T^3)/n
  % then the rate of entropy production is
  % = (tr(P)*tr(Pinv)-9)/(1.5*tau_0*sqrt(m_s*T^3))
  [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
    = get_primVars10(c, state_coef,species);
  m_s = get_species_mass(c, species);
  trP = P11 + P22 + P33;
  PA11 = P22.*P33 - P23.*P23;
  PA22 = P11.*P33 - P13.*P13;
  PA33 = P11.*P22 - P12.*P12;
  PA12 = P23.*P13 - P12.*P33;
  PA13 = P12.*P23 - P22.*P13;
  detP = PA11.*P11 + PA12.*P12 + PA13.*P13;
  trPinv = (PA11 + PA22 + PA33)./detP;
  n = (1./m_s)*rho;
  [tau_0] = get_base_iso_period(c, species);
  switch c.s.params.iso_period_type
  case 'constant'
    disp(['using constant tau_0 = ' num2str(tau_0)]);
    Sst = (trP.*trPinv-9).*n.*(2./(3*tau_0));
  case 'det'
    disp(['using detT with tau_0 = ' num2str(tau_0)]);
    detT = detP./(n.^3);
    Sst = (2./(3*tau_0*sqrt(m_s)))*(trP.*trPinv-9.)./sqrt(detT);
  case 'trace'
    disp(['using trT with tau_0 = ' num2str(tau_0)]);
    T = (m_s/3.)*(trP./rho);
    Sst = (2./(3*tau_0*sqrt(m_s)))*(trP.*trPinv-9.)./(T.^(3./2.));
  otherwise
    error(['unsupported iso_period_type: ' species]);
  end
end

function val = sample_Pi11(c, state_coef); val = get_Psij(c, state_coef, 'i',1,1); end
function val = sample_Pi22(c, state_coef); val = get_Psij(c, state_coef, 'i',2,2); end
function val = sample_Pi33(c, state_coef); val = get_Psij(c, state_coef, 'i',3,3); end
function val = sample_Pi12(c, state_coef); val = get_Psij(c, state_coef, 'i',1,2); end
function val = sample_Pi13(c, state_coef); val = get_Psij(c, state_coef, 'i',1,3); end
function val = sample_Pi23(c, state_coef); val = get_Psij(c, state_coef, 'i',2,3); end
function val = sample_Pe11(c, state_coef); val = get_Psij(c, state_coef, 'e',1,1); end
function val = sample_Pe22(c, state_coef); val = get_Psij(c, state_coef, 'e',2,2); end
function val = sample_Pe33(c, state_coef); val = get_Psij(c, state_coef, 'e',3,3); end
function val = sample_Pe12(c, state_coef); val = get_Psij(c, state_coef, 'e',1,2); end
function val = sample_Pe13(c, state_coef); val = get_Psij(c, state_coef, 'e',1,3); end
function val = sample_Pe23(c, state_coef); val = get_Psij(c, state_coef, 'e',2,3); end

function detP=sample_detPi(c, state_coef); detP = sample_detPs(c, state_coef, 'i'); end
function detP=sample_detPe(c, state_coef); detP = sample_detPs(c, state_coef, 'e'); end
function detT=sample_detTi(c, state_coef); detT = sample_detTs(c, state_coef, 'i'); end
function detT=sample_detTe(c, state_coef); detT = sample_detTs(c, state_coef, 'e'); end

function [detP, p, rho] = sample_detPs(c, state_coef, species);
  species_model=get_species_model(c, species);
  switch species_model
  case 10
    [rho, u1, u2, u3, P11, P12, P13, P22, P23, P33] ...
       = get_primVars10(c, state_coef, species);
    detP = get_determinant(P11, P12, P13, P22, P23, P33);
    p = (P11 + P22 + P33)*(1./3.);
  case 5
    [rho, u1, u2, u3, p] = get_primVars05(c, state_coef, species);
    detP = p.^3;
  otherwise
    error(['unsupported model: ' species_model]);
  end
end

function [detT, p, n] = sample_detTs(c, state_coef, species);
  [detP, p, rho] = sample_detPs(c, state_coef, species);
  species_model=get_species_model(c, species);
  ms = get_species_mass(c, species);
  n = rho/ms;
  detT = detP./(n.^3);
end

% handle mirroring and wrap-around without actually padding the array.
function out = dx_state(c, state, space_order,dx,flip)
  assert(space_order >= 1);
  assert(space_order <= 3);
  dims=size(state);
  assert(dims(3)==get_kmax(space_order));
  out = zeros(dims);
  paddims = dims;
  paddims(1) = 1;
  kmax=size(state,3);

  % determine ghost cell contents
  % based on enforced symmetries and the periodic boundary conditions
  %
  global enforced_symmetry;
  leftpad = zeros(paddims);
  rghtpad = zeros(paddims);
  enforced_x_axis_symmetry = bitand(c.s.params.enforced_symmetry,1);
  if(enforced_x_axis_symmetry)
    reversed_components = [1, -1,1, -1,1,1]; % 1st, 2nd, 3rd order components
    left_sgns = reversed_components*flip(1);
    rght_sgns = reversed_components*flip(2);
    for i=1:kmax
      leftpad(1,:,i,:) = left_sgns(i)*state(  1,:,i,:);
      rghtpad(1,:,i,:) = rght_sgns(i)*state(end,:,i,:);
    end
  else % periodic
    % for the periodic case simply copy wrapped array elements
    leftpad(1,:,:,:) = state(end,:,:,:);
    rghtpad(1,:,:,:) = state(  1,:,:,:);
  end

  % take centered differences of all components
  out(2:end-1,:,:,:) = (state(3:end,:,:,:)-state(1:end-2,:,:,:))/(2.*dx);
  out(1,:,:,:) = (state(2,:,:,:)-leftpad)/(2.*dx);
  out(end,:,:,:) = (rghtpad - state(end-1,:,:,:))/(2.*dx);
  if(space_order>=3)
    % correct first component for third-order accuracy
    out(2:end-1,:,1,:) = out(2:end-1,:,1,:) ...
      - (state(3:end,:,5,:)-state(1:end-2,:,5,:))*(sqrt(5)/dx);
    out(1,:,1,:) = out(1,:,1,:) ...
      - (state(2,:,5,:)-leftpad(1,:,5,:))*(sqrt(5)/dx);
    out(end,:,1,:) = out(end,:,1,:) ...
      - (rghtpad(1,:,5,:) - state(end-1,:,5,:))*(sqrt(5)/dx);
  end
end

function out = dy_state(c, state, space_order,dy,flip)
  assert(all(flip.*flip==[1 1]));
  assert(space_order >= 1);
  assert(space_order <= 3);
  dims=size(state);
  out = zeros(dims);
  paddims = dims;
  paddims(2) = 1;
  kmax=size(state,3);
  %
  % determine ghost cell sources and negations
  % based on boundary conditions and enforced symmetries
  CONDUCTING_WALL=1; PERIODIC=2;
  BCs = CONDUCTING_WALL;
  enforced_y_axis_symmetry = bitand(c.s.params.enforced_symmetry,2);
  % determine ghost cell contents
  % based on enforced symmetries and the periodic boundary conditions
  low_pad = zeros(paddims);
  highpad = zeros(paddims);
  %
  reversed_components = [1, 1,-1, -1,1,1];
  low__sgns = reversed_components*flip(1);
  high_sgns = reversed_components*flip(2);
  %
  if(enforced_y_axis_symmetry)
    % low end obeys mirror symmetry
    %disp('enforcing vertical mirror symmetry');
    for i=1:kmax
      low_pad(:,1,i,:) = low__sgns(i)*state(:,  1,i,:);
      highpad(:,1,i,:) = high_sgns(i)*state(:,end,i,:);
    end
  else
    if(BCs==PERIODIC)
      % for the periodic case simply copy wrapped array elements
      low_pad(:,1,:,:) = state(:,end,:,:);
      highpad(:,1,:,:) = state(:,  1,:,:);
    else % CONDUCTING WALL
      for i=1:kmax
        %low_pad(:,1,i,:) = low__sgns(i)*state(:,  1,i,:);
        %highpad(:,1,i,:) = high_sgns(i)*state(:,end,i,:);
        low_pad(:,1,i,:) = low__sgns(i)*state(:,  2,i,:);
        highpad(:,1,i,:) = high_sgns(i)*state(:,end-1,i,:);
      end
    end
  end

  % take centered differences of all components
  %
  out(:,2:end-1,:,:) = (state(:,3:end,:,:)-state(:,1:end-2,:,:))/(2.*dy);
  out(:,1,:,:) = (state(:,2,:,:)-low_pad)/(2.*dy);
  out(:,end,:,:) = (highpad - state(:,end-1,:,:))/(2.*dy);
  if(space_order>=3)
    % correct first component for third-order accuracy
    out(:,2:end-1,1,:) = out(:,2:end-1,1,:) ...
      - (state(:,3:end,6,:)-state(:,1:end-2,6,:))*(sqrt(5)/dy);
    out(:,1,1,:) = out(:,1,1,:) ...
      - (state(:,2,6,:)-low_pad(:,1,6,:))*(sqrt(5)/dy);
    out(:,end,1,:) = out(:,end,1,:) ...
      - (highpad(:,1,6,:) - state(:,end-1,6,:))*(sqrt(5)/dy);
  end
end

function qm = get_qm_ratio(c, species)
  switch species
  case 'i'
    qm = 1./c.s.params.ion_mass;
  case 'e'
    qm = -1./c.s.params.elc_mass;
  otherwise
    error(['unsupported species: ' species]);
  end
end

function dtus3 = sample_dtus3(c, state_coef, species)
  E3 = sample_E3(c, state_coef);
  [Bxus3, rho_s] = sample_Bxus3(c, state_coef, species);
  qm = get_qm_ratio(c, species);
  dtus3 = qm.*(E3-Bxus3) - sample_divPs3(c,state_coef,species)./rho_s;
end

function dtui3 = sample_dtui3(c, state_coef)
  dtui3 = sample_dtus3(c, state_coef, 'i');
end

function dtue3 = sample_dtue3(c, state_coef)
  dtue3 = sample_dtus3(c, state_coef, 'e');
end

function data = sample_dtJ(c, state_coef);
  %global space_order;
  %global rho_i Ni B E;
  %global ion_mass;
  %global dx dy;

  % Suffixes indicate one of three ways of representing states:
  %
  % state_coef = sc = lc: legendre coefficients
  % state_vals = sv = qv: quadrature values.
  % plot_vals  = pv = rg: regular grid samples
  %
  % high-order numerical differentation needs the lc representation:
  %   div_state_lc = div_stateTensor(state_lc,space_order),
  % arithmetic operation needs the qv or rg representation, and
  % for output we use the rg representation.
  %
  % conversions supported are qv <-> lc -> rg, i.e.,:
  %   qv -> lc
  %   lc -> qv
  %   lc -> rg
  %
  % Essentially we need to stick with the lc/qv representation
  % until we are done with differentiation, at which point we
  % can move to rg representation.  So if we want to calculate
  % the *curl* of the Ohm's law terms, we will have to delay
  % moving to the rg representation.
  %
  % conversion works by:
  %
  % state_lc = project_onto_basis_cart2(state_qv,space_order);
  % state_qv = sample_state2(state_lc, space_order, 2);
  % state_rv = sample_state2(state_lc, space_order, 1);

  stateIndices_s = c.s.stateVarMap.stateIndices_s;
  how_sample = c.s.sample_rate;
  % the multiplier
  %
  rho_i_lc = state(:,:,:,stateIndices_s.rho_i);
  rho_i_qv = sample_state_cart2(rho_i_lc,how_sample);
  rho_qv = 2.0*rho_i_qv;
  mass_prod = ion_mass*ion_mass;
  mass_prod_ovr_rho = mass_prod./rho_qv;

  % for the source term
  %
  E_lc = state(:,:,:,stateIndices_s.E);
  E_qv = sample_state_cart2(E_lc,how_sample);
  B_lc = state(:,:,:,stateIndices_s.B);
  B_qv = sample_state_cart2(B_lc,how_sample);

  % ions
  %
  % flux term
  %
  Ni_state=state(:,:,:,stateIndices_s.Ni([3 5]));
  Mi3_flux_lc = div_stateTensor2(Ni_state,space_order,dx,dy);
  Mi3_flux_qv = sample_state_cart2(Mi3_flux_lc, how_sample);

  J3_t_qv = 2.0*Mi3_t_qv/ion_mass;

  % DivJflux_lc
  %
  %Mi3_state = state(:,:,:,stateIndices_s.Mi(3));
  %rho_i_state = state(:,:,:,stateIndices_s.rho_i);
  %u_qv = sample_u_val(Mi_state,rho_i_state,space_order,2);
  %J3_state = (2.0/ion_mass)*state(:,:,:,stateIndices_s.Mi(3));
  %J3_qv = sample_state2(J3_state,space_order,2);
  %uJ_plus_Ju_qv = twiceSymmetricTensorProduct(u_qv, J3_qv);
  %Jflux_lc = project_onto_basis_cart2(uJ_plus_Ju_qv,space_order);
  %DivJflux_lc = div_stateTensor2(Jflux_lc,space_order,dx,dy);
  %DivJflux_qv = sample_state2(DivJflux_lc, space_order, 2);

  J_t_term_qv = tupleTimesScalar(J3_t_qv,mass_prod_ovr_rho);
  % why do I fail to get agreement when I include DivJflux_qv?
  intertial_term_qv = J_t_term_qv; % + DivJflux_qv;

  % this step would be unnecessary if we had worked with
  % rg values, but I wanted to retain the possibility of converting
  % the results to lc representation so we could take the curl
  data = convert_qv_to_rv(intertial_term_qv, space_order);
end


