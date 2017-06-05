
function [mx,my,xlow,xhigh,ylow,yhigh,method1,meqn,nout,mx_out,my_out,datafmt]=get_plot_params2(outputdir);

  [nout,tfinal,dtv,cflv,nv,...
     method,meqn,mx,my,mbc,xlow,xhigh,ylow,yhigh,...
     mrestart,nstart,datafmt]...
   = read_dogpack_parameters2([outputdir '/dogpack.data']);

  method1=method(1);
  mx_out = mx*method1;
  my_out = my*method1;

end
