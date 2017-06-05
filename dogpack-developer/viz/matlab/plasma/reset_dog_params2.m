
function [xlow, xhigh, ylow, yhigh, mx, my] ...
= reset_dog_params2(mx, my, domain_scaling, ...
    enforced_symmetry, sym_comp_dims, use_centered_cells);

    XBIT=1;
    YBIT=2;

    % handle extra bits in enforced_symmetry for
    % backwards compatibility

    % ...

    xhigh = 4.*pi*domain_scaling;
    yhigh = 2.*pi*domain_scaling;
    xlow = -xhigh;
    ylow = -yhigh;

    if(bitand(enforced_symmetry, XBIT))
        if(bitand(sym_comp_dims, XBIT))
            if(bitand(use_centered_cells, XBIT))
                dx = 2.*xhigh/(2*mx-1);
                xlow = -dx/2.;
            else
                xlow = 0;
            end
        end
        else if(bitand(mx,1)==0)
        % mx is even
            xlow = 0.;
            mx = mx/2;
        else
        % mx is odd
            dx = 2*xhigh/mx;
            xlow = -dx/2.;
            mx = (mx+1)/2;
        end
    end

    if(bitand(enforced_symmetry, YBIT))
        if(bitand(sym_comp_dims, YBIT))
            if(bitand(use_centered_cells, YBIT))
                dy = 2.*yhigh/(2*my-1);
                ylow = -dy/2.;
            else
                ylow = 0;
            end
        end
        else if(bitand(my,1)==0)
        % my is even
            ylow = 0.;
            my = my/2;
        else
        % my is odd
            dy = 2*yhigh/my;
            ylow = -dy/2.;
            my = (my+1)/2;
        end
    end
end
