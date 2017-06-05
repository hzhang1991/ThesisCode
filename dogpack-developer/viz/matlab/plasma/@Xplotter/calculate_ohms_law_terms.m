
function [mu3i, inertial_term, divPs_term, Ps13_x_term, Ps23_y_term, ...
    resistive_term, residual] = ...
  calculate_ohms_law_terms(xplotter,t, rho_i, M3i,  P13i_x, P23i_y, E3, ...
    species, spc_mass, params)

  % calculate ion momentum (ohm's law) quantities
  charge = 1;
  if(species=='e')
    charge = -1;
  end
  spc_charge_to_mass_ratio = charge/spc_mass;
  ns_q = rho_i.*spc_charge_to_mass_ratio;
  %spc_mass
  Ps13_x_term = P13i_x./ns_q;
  Ps23_y_term = P23i_y./ns_q;
  % n_i = rho_i./spc_mass;
  mu3i = M3i./ns_q;

  inertial_term = calculate_inertial_term(xplotter, t, mu3i);
  % shouldn't I be dividing the resistive term by n_i?
  resistive_term = calculate_resistive_term(mu3i,params);

  % calculate Ohm's law terms
  divPs_term = Ps13_x_term + Ps23_y_term;
  % probably there is a better way to do this
  %if(species=='e')
  %  P13i_x_term =  -P13i_x_term;
  %  P23i_y_term =  -P23i_y_term;
  %  divPs_term=-divPs_term;
  %  inertial_term=-inertial_term;
  %  resistive_term=-resistive_term;
  %end

  residual = E3-inertial_term-divPs_term-resistive_term;
end

function resistive_term = calculate_resistive_term(mu3i,params)
  resistive_term = zeros(size(mu3i));
  if(params.slowing_rate_slope > 0)
  % anomalous resistivity
    abs_mu3i = abs(mu3i);
    ind = find(abs_mu3i>params.trigger_mui3);
    local_slowing_rate = zeros(size(mu3i));
    local_slowing_rate(ind) = params.slowing_rate_slope...
       *(abs_mu3i(ind)-params.trigger_mui3);
    resistive_term = local_slowing_rate.*mu3i;
  elseif(params.slowing_period > 0)
  % case of constant resistivity
    slowing_rate = 1./params.slowing_period;
    resistive_term = mu3i*slowing_rate;
  end
end

