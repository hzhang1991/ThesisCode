% flip = 1 for vector, -1 for pseudo-vector
function vec_norm = plot_vector(xc_in,yc_in,vecx_in,vecy_in,scale,...
   flip, enforced_symmetry);

    if(~exist('scale','var'))
       scale=.4;
    end
    % set skip so that
    xnum = size(xc_in,1);
    xnum_desired = 20;
    xskip_f = xnum/xnum_desired;
    xskip = ceil(xskip_f);
    xstart = ceil(xskip_f/2.);

    xc=xc_in(xstart:xskip:end,xstart:xskip:end);
    yc=yc_in(xstart:xskip:end,xstart:xskip:end);
    vecx=vecx_in(xstart:xskip:end,xstart:xskip:end);
    vecy=vecy_in(xstart:xskip:end,xstart:xskip:end);

    % automatic vector scaling doesn't seem to work according to
    % the documentation when scale~=0, so I'll do it myself.
    %
    if(0)
      % compute L2 norm (rms) of vector
      vec2 = (vecx.^2 + vecy.^2);
      vec_norm = sqrt(mean(mean(vec2)));
    else
      % scale vectors using the max norm
      % compute the infinity norm (max norm)
      vec_norm = max(max(max(abs(vecx))),max(max(abs(vecy))));
      disp(['infinity norm of 2D vector field: ' num2str(vec_norm)]);
    end
    %
    %vec_scale=scale/(1.+vec_norm);
    vec_scale = scale/vec_norm;
    if(scale~=0) %if(scale~=0)
      % rescale the magnitudes by rms magnitude for display
      vecx = vecx.*vec_scale;
      vecy = vecy.*vec_scale;
    end

    % center arrows rather than placing tail at position
    xc=xc-.5*vecx;
    yc=yc-.5*vecy;

    if(bitand(enforced_symmetry,3)==0)
      cx=xc;
      cy=yc;
      vx=vecx;
      vy=vecy;
      %h=quiver(xc,yc,vecx,vecy,0);
    else
      if(bitand(enforced_symmetry,3)==3)
        if(flip==1)
          vx = [-vecx(end:-1:1,end:-1:1), vecx(:,end:-1:1); ...
                -vecx(end:-1:1,:),vecx];
          vy = [-vecy(end:-1:1,end:-1:1),-vecy(:,end:-1:1); ...
                 vecy(end:-1:1,:),vecy];
        else
          vx = [-vecx(end:-1:1,end:-1:1),-vecx(:,end:-1:1); ...
                 vecx(end:-1:1,:),vecx];
          vy = [-vecy(end:-1:1,end:-1:1), vecy(:,end:-1:1); ...
                -vecy(end:-1:1,:),vecy];
        end
        cx = [-xc(end:-1:1,end:-1:1), xc(:,end:-1:1); ...
              -xc(end:-1:1,:),xc];
        cy = [-yc(end:-1:1,end:-1:1),-yc(:,end:-1:1); ...
               yc(end:-1:1,:),yc];
      elseif(bitand(enforced_symmetry,3)==2)
        if(flip==1)
          vx = [ vecx(:,end:-1:1); vecx];
          vy = [-vecy(:,end:-1:1); vecy];
        else
          vx = [-vecx(:,end:-1:1); vecx];
          vy = [ vecy(:,end:-1:1); vecy];
        end
        cx = [ xc(:,end:-1:1); xc];
        cy = [-yc(:,end:-1:1); yc];
      elseif(bitand(enforced_symmetry,3)==1)
        error('unimplemented');
      else
        error('impossible');
      end
    end
    h=quiver(cx,cy,vx,vy,0);
    set(h,'Color','black');

end

