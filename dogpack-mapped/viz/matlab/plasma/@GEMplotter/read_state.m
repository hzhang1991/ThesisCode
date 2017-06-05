
function [state,time] = read_state(c,frame_number,components)
  params = c.s.params;
  if(~exist('components','var'))
    components=1:params.meqn;
  end
  [state,time]=read_state2_cart(params.datafmt, params.outputdir, ...
    frame_number, 'q',...
    params.plot_mx, params.plot_my, ...
    params.meqn, get_kmax(params.space_order), components);
end

