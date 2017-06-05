%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                  mx, my:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh]:  min/max values of grid
%%                    meqn:  number of equations
%%                    maux:  number of aux components
%%                   meth1:  spatial order of accuracy
%%
%%   Grid information:
%%       (xc,yc): grid points (cell centers), size = (mx,my)
%%       (xl,yl): grid points (lower left cell corners), size = (mx+1,my+1)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,maux)
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,meqn)
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear qex;
% for j=1:my
%   for i=1:mx

%     qex(i,j) = sin( 2*pi*(xc(i,j)-pi*time) )* exp(-4.0*yc(i,j)^2);

%       qex(i,j) = 2.0*cos(2.0*pi*xc(i,j) - yc(i,j)*time )*sin(pi*yc(i,j) );
%       rr = sqrt((xc(i,j)-0.4)^2 + (yc(i,j)-0.5)^2);
%       xold = mod( xc(i,j) - yc(i,j)*time, 1 );
%       rr = sqrt((xold-0.4)^2 + (yc(i,j)-0.5)^2);
%       if (rr<0.3)
%         qex(i,j) = cos(5.0/3.0*pi*rr).^6;
%       else
%         qex(i,j) = 0.0;
%       end
%   end
% end

  clear qex;
  for j=1:my
    for i=1:mx
      rr = sqrt((xc(i,j)-0.4)^2 + (yc(i,j)-0.5)^2);
      if (rr<0.3)
        qex(i,j) = (cos(5.0/3.0*pi*rr)).^6;
      else
        qex(i,j) = 0.0;
      end
    end
  end


  if (abs(time)<=1.0e-12 || abs(time-1)<=1.0e-12 || abs(time-10) <= 1.0e-12)
    err = reshape(abs(qex-qsoln(:,:,1)),mx*my,1);
    err_scale = reshape(abs(qex),mx*my,1);
  
    err_rel = norm(err,2)/norm(err_scale,2);

    err_rel = norm(qsoln-qex,2)/norm(qex,2);

    disp([num2str(err_rel,'%0.10e')]);
  end
