
function set_axes(domain_scaling,xlow,xhigh,ylow,yhigh)
    use_defaults=0;
    if(nargin<5)
      use_defaults=1;
      xlow=-4*pi;
      xhigh=4*pi;
      ylow=-2*pi;
      yhigh=2*pi;
    end
    xlow = xlow * domain_scaling;
    ylow = ylow * domain_scaling;
    xhigh = xhigh * domain_scaling;
    yhigh = yhigh * domain_scaling;

    axis on; box on; grid off;
    axis('equal');
    axis([xlow xhigh ylow yhigh]);
    xtick = (-12:4:12)*domain_scaling;
    ytick = (-6:2:6)*domain_scaling;
    if(use_defaults)
       set(gca,'xtick',xtick);
       set(gca,'ytick',ytick);
    end
    axes_left=0.05;
    axes_bottom=0.08;
    axes_width = .9;
    axes_height = .82;
    set(gca,...
      'fontsize',24,...
      'fontweight','bold',...
      'linewidth',2,...
      'position',[axes_left axes_bottom axes_width axes_height]);
end

