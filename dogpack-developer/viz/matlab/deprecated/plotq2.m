% generic rudimentary routine to plot output data
%
function plotq2(nframe, outputdir_in)
  global outputdir
  if(nargin>=2)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end
  
  format long e;

  q = readq2(nframe, outputdir)

  %[meqn,mx,my,method1,nout,xl,yl,xc,yc,datafmt]=readgrid(outputdir);
  [mx,my,xlow,xhigh,ylow,yhigh,method1,meqn,nout,mx_out,my_out,datafmt]...
    = get_plot_params2(outputdir);
  [dx,dy,xl,yl,xc,yc]=get_grid2(mx,my,xlow,xhigh,ylow,yhigh,method1);
  %
  axis_array=[xlow xhigh ylow yhigh];

    
  assert(numel(nframe)==1,['frame number is not a scalar: ' num2str(nframe)]);

      [qfull,time]=readq2(n1,outputdir,datafmt,meqn,mx,my,method1,mc);

      for i=1:length(mc)
         % perhaps the user should be able to override this.
         plot_scalar2(...
            mc(i),qfull(:,:,i),xl,yl,time,axis_array,['q(' num2str(mc(i))]);
         if(isa(after_plot_scalar, 'function_handle'))
            % This function call is a generic message,
            % so broadcast all relevant information.
            after_plot_scalar;
         end
      end
      if(isa(after_plots, 'function_handle'))
         after_plots;
      end
    end
  end
  disp(' ')
end

function plot_scalar2(fig,qvals,xl,yl,time,axis_array,titleStr);

    if(ishandle(fig))
      set(0,'CurrentFigure',fig);
    else
     figure(fig);
    end
    clf;
    % to accomodate pcolor's idiosyncrasies...
    qplot=zeros(size(qvals)+1);
    qplot(1:end-1,1:end-1)=qvals;
    pcolor(xl,yl,qplot);
    shading flat;
    colormap('jet');
    hold on;
    axis on; box on; grid off;
    axis('equal');
    axis(axis_array);
    % could create a hook like "after_plot_scalar"
    %axis([0 4*pi 0 2*pi]);
    %set(gca,'xtick',-12:4:12);
    %set(gca,'ytick',-6:2:6);
    set(gca,'fontsize',16);
    t1 = title([titleStr ' at t = ',num2str(time) ' ']); 
    set(t1,'fontsize',16);
    %caxis([0,1]);
    c1=colorbar;
    set(c1,'fontsize',16);
end
