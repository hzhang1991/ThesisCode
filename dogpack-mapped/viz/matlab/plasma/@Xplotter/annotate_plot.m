
function annotate_plot(xplotter,params,plot_struct,fonts, ...
  using_resistive_term, upper_xshift)

  annotate_plot_PoP(xplotter,params,plot_struct,fonts, ...
  using_resistive_term, upper_xshift)
end

function annotate_plot_PoP(xplotter,params,plot_struct,fonts, ...
  using_resistive_term, upper_xshift)

  limits = axis();
  xpos = limits(1)+.02*(limits(2)-limits(1));
  ypos = limits(3)+.16*(limits(4)-limits(3));
  %text(xpos, ypos, ['(rescaled problem) '],'fontsize',fonts.annotationfontsize,'fontweight','bold');
  %text(xpos, ypos, ['(GEM-like scaling) '],'fontsize',fonts.annotationfontsize,'fontweight','bold');

  model_string = get_model_string(Xplotter(),params, plot_struct.do_rescale_times);
  xpos = limits(1)+.02*(limits(2)-limits(1));
  ypos = limits(3)+.05*(limits(4)-limits(3));
  text(xpos, ypos, [model_string ' '],'fontsize',fonts.annotationfontsize,'fontweight','bold');
end

function annotate_plot_full(xplotter,params,plot_struct,fonts, ...
  using_resistive_term, upper_xshift)

  method_string = get_method_string(Xplotter(),params);
  %method_string = plot_struct.outputdir_shouldname;
  model_string = get_model_string(Xplotter(),params, plot_struct.do_rescale_times);
  full_problem_string = get_full_problem_string(Xplotter(),params,...
    plot_struct.do_rescale_lengths);
  %model_problem_string = get_model_problem_string(Xplotter(),params,...
  %  plot_struct.do_rescale_lengths);
  limits = axis();
  xpos = limits(1)+.02*(limits(2)-limits(1));
  ypos = limits(3)+.09*(limits(4)-limits(3));
  text(xpos, ypos, [model_string ' '],'fontsize',fonts.annotationfontsize,'fontweight','bold');
  xpos = limits(1)+.02*(limits(2)-limits(1));
  ypos = limits(3)+.04*(limits(4)-limits(3));
  text(xpos, ypos, full_problem_string,'fontsize',fonts.annotationfontsize,'fontweight','bold');
  xpos = limits(1)+upper_xshift*(limits(2)-limits(1));
  ypos = limits(3)+.94*(limits(4)-limits(3));
  text(xpos, ypos, method_string,'fontsize',fonts.annotationfontsize,'fontweight','bold');
  ypos2 = limits(3)+.90*(limits(4)-limits(3));
  if(using_resistive_term)
    slowing_string = get_slowing_string(params, ...
      plot_struct.do_rescale_times,plot_struct.do_rescale_fluxes);
    text(xpos, ypos2, slowing_string,'fontsize',fonts.annotationfontsize,'fontweight','bold');
  end
  % say if it blew up
  blewup = test_blewup(Xplotter(), params.outputdir);
  if(blewup)
    xpos = limits(1)+.85*(limits(2)-limits(1));
    ypos = limits(3)+.40*(limits(4)-limits(3));
    text(xpos, ypos, '(blew up)','fontsize',fonts.annotationfontsize,'fontweight','bold');
  end
end

