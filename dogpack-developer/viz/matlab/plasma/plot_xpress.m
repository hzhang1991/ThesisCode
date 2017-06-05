function plot_xpress(outputdir_in, final_time, figure_shift)
  setarg_outputdir
  if(~exist('figure_shift')) figure_shift=20; end
  format long e
  
  %%%%%%%%%%%%%%%%%%%%
  %%% read in the data
  %%%%%%%%%%%%%%%%%%%%

  fid = fopen([outputdir '/xpoint_gas_i.dat']);
  xpoint_gas = fscanf(fid,'%e',[5 inf]);
  status = fclose(fid);
  %
  t = xpoint_gas(1,:)';
  rho_i = xpoint_gas(2,:)';
  M3i = xpoint_gas(3,:)';
  P13i_x = xpoint_gas(4,:)';
  P23i_y = xpoint_gas(5,:)';
  nt = length(t);
  
  fid = fopen([outputdir '/xpoint_E3.dat']);
  xpoint_E3 = fscanf(fid,'%e',[2 inf]);
  status = fclose(fid);
  %
  E3 = xpoint_E3(2,:)';

  %%%%%%%%%%%%%%%%%%%%
  %%% calculate quantities
  %%%%%%%%%%%%%%%%%%%%
  
  % calculate ion momentum (ohm's law) quantities
  divPi = P13i_x+P23i_y;
  u3i = M3i./rho_i;
  u3i_t_middle = (u3i(3:end)-u3i(1:end-2)) ...
              ./ (  t(3:end)-  t(1:end-2));
  % extrapolate for boundary values
  u3i_t_left = 2*u3i_t_middle(1) - u3i_t_middle(2);
  u3i_t_rght = 2*u3i_t_middle(end) - u3i_t_middle(end-1);
  u3i_t = zeros(size(u3i));
  u3i_t(1) = u3i_t_left;
  u3i_t(end) = u3i_t_rght;
  u3i_t(2:end-1) = u3i_t_middle;
  
  % get ion_mass
  parameters_ini = [outputdir '/parameters.ini'];
  params = read_params_from_ini(parameters_ini, ...
    'mass_ratio', 'mx', 'my', 'model_name', ...
    'ion_iso_period', 'elc_iso_period');
  mass_ratio = params.mass_ratio;
  elc_mass = 1./(mass_ratio+1.);
  ion_mass = mass_ratio*elc_mass;
  %ion_mass = 0.5; % should calculate from parameters.
  
  % calculate Ohm's law terms
  n_i = rho_i/ion_mass;
  divPi_term    = divPi./n_i    ;
  inertial_term = ion_mass*u3i_t;
  resistive_term = E3-inertial_term-divPi_term;

  if(accumulate)
   int_E3             = integrate(t,E3);
   int_divPi_term     = integrate(t,divPi_term);
   int_inertial_term  = integrate(t,inertial_term);
   int_resistive_term = integrate(t,resistive_term);
  end
  
  effective_resistivity = resistive_term./u3i;

  %%%%%%%%%%%%%%%%%%%%
  %%% plot data
  %%%%%%%%%%%%%%%%%%%%
  if(~exist('final_time'))
    final_time = t(end);
  end
  
% % plot Ohm's law terms at X-point
% %
% figure(1+figure_shift);
% clf;
% hold on;
% p1=plot(t,-E3             ,'b-'); set(p1,'linewidth',2);
% p1=plot(t,-divPi_term     ,'g-'); set(p1,'linewidth',2);
% p1=plot(t,-inertial_term  ,'k-'); set(p1,'linewidth',3);
% p1=plot(t,-resistive_term ,'m-'); set(p1,'linewidth',2);
% hold off;
% set(gca,'fontsize',16);
% t1=title([ 'Ohm''s law terms at the X-point [DoGPack]']);
% set(t1,'fontsize',32);
% axis([t(1) final_time -.3 .3]);
% axis on; box on; grid off;
% hold off
% %set(gca,'xtick',0:100:500);
% %set(gca,'ytick',0:0.01:0.15);
% set(gca,'plotboxaspectratio',[1.5 1 1]);
% l1=legend(...
%    'E_3', ...
%    'div(P_i)_3/(e n_i)',...
%    '(m_i/e)d_t (u_i)_3',...
%    'resistive term'...
%     );
% set(l1,'location','northeast');
% set(l1,'fontsize',26);
  
 % plot cumulation integrals
 %
 figure(2+figure_shift);
 clf;
 hold on;
 p1=plot(t,-int_E3             ,'b-' ); set(p1,'linewidth',3);
 p1=plot(t,-int_divPi_term     ,'r--'); set(p1,'linewidth',2);
 p1=plot(t,-(ion_mass*u3i)     ,'k:' ); set(p1,'linewidth',2);
 p1=plot(t,-int_inertial_term  ,'k-.' ); set(p1,'linewidth',2);
 p1=plot(t,-int_resistive_term ,'m'); set(p1,'linewidth',2);
 hold off;
 set(gca,'fontsize',16);
 mesh_string = ['(' num2str(params.mx) 'x' num2str(params.my) ' mesh)' ];
 t1=title([ 'accumulation integral of Ohm''s law terms at the X-point ' ]);
 set(t1,'fontsize',24);
 axis on; box on; grid off;
 xlim([t(1), final_time]);
 ylim([-.5 2.5]);
 hold off
 %set(gca,'xtick',0:100:500);
 %set(gca,'ytick',0:0.01:0.15);
 set(gca,'plotboxaspectratio',[1.5 1 1]);
 limits = axis();
 if(strcmp(params.model_name,'p05'))
   model_string = '5-moment pair plasma'
 elseif(strcmp(params.model_name,'p10'))
   ion_iso_period = params.ion_iso_period;
   assert(params.ion_iso_period==params.elc_iso_period);
   if(ion_iso_period<0)
     model_string = ['10-moment pair plasma, no isotropization'];
   elseif(ion_iso_period==0)
     model_string = ['10-moment pair plasma, instantaneous isotropization'];
   else
     model_string = ['10-moment pair plasma, isotropization period = ' ...
                   num2str(ion_iso_period)];
   end
 else
   error(['invalid model: ' params.model_name]);
 end
 xpos = limits(1)+.02*(limits(2)-limits(1));
 ypos = limits(3)+.12*(limits(4)-limits(3));
 text(xpos, ypos, model_string,'fontsize',20);
 xpos = limits(1)+.02*(limits(2)-limits(1));
 ypos = limits(3)+.05*(limits(4)-limits(3));
 text(xpos, ypos, mesh_string,'fontsize',22);
 %l1=legend(...
 %   '$$\int_0^t E_3$$', ...
 %   '$$\int_0^t (\nabla. P_i)_3/(e n_i)$$',...
 %   '$$(m_i/e) (u_i)_3$$',...
 %   '$$\int_0^t (m_i/e)\partial_t (u_i)_3$$',...
 %   '$$\int_0^t \hbox{resistive term}$$'...
 %    );
 %set(l1,'Interpreter','latex');
 l1=legend(...
    '$\int_0^t E_3$', ...
    '$\int_0^t \frac{(\nabla. P_i)_3}{e n_i}$',...
    '$\frac{m_i}{e} (u_i)_3$',...
    '$\int_0^t \frac{m_i}{e}\partial_t (u_i)_3$',...
    '$\int_0^t \hbox{residual}\ \ $'...
     );
 set(l1,'Interpreter','latex');
 %text('Interpreter','latex',...
 % 'String','$$\int_0^x\!\int_y dF(u,v)$$',...
 %  'Position',[.5 .5],...
 %   'FontSize',16);
 set(l1,'location','northwest');
 set(l1,'fontsize',24);
 xlabel('time','fontsize',20)
  
% % plot effective resistivity
% %
% figure(3+figure_shift);
% clf;
% hold on;
% p1=plot(t,effective_resistivity ,'k-'); set(p1,'linewidth',2);
% hold off;
% set(gca,'fontsize',16);
% t1=title(['effective "resistivity" (resistive term)/(u_i)_3']);
% set(t1,'fontsize',24);
% axis([t(1) final_time -.1 .3]);
% % axis 'auto y';
% axis on; box on; grid off;
% hold off
% %set(gca,'xtick',0:100:500);
% %set(gca,'ytick',0:0.01:0.15);
% set(gca,'plotboxaspectratio',[1.5 1 1]);
  
end
  
function out=integrate(t,f)
  fh = (f(2:end) + f(1:end-1))*0.5;
  dt = t(2:end)-t(1:end-1);
  out=zeros(size(f));
  out(2:end) = cumsum(fh.*dt);
end
