
function method_string = get_method_string(xplotter,params)
  method_string = ...
    ['computational grid: ' num2str(params.mx) 'x' num2str(params.my)];
  assert(params.time_order==params.space_order);
  order = params.time_order;
  max_cfl = 1./(2*order-1.);
  method_string = [method_string ...
    ', order = ' num2str(params.space_order) ...
    ... % ', cfl = '   num2str(params.cflv(2)) ... % comment out?
    ... % ' (max cfl = ' num2str(max_cfl) ')' ... % comment out?
    ', cfl/cfl_{max} = '   sprintf('%.2f',params.cflv(2)/max_cfl) ... % comment out?
    ' ' ];
    %', cfl = '   num2str(params.cflv(2)) '/' num2str(max_cfl) ...
end

