%MAIN_DRIVER2D.
%
% This script can be used to run a refinement study for a solution computed by
% DoGPack.
%
% This script has the ability to plot pointwise errors and postprocessed 
% pointwise errors.
%
% The variables you have at your disposal include:
%
%       runname : a name used to identify what run you are looking at.
%                 Data files should be found in 'output[runname]000/,
%                 output[runname]001/, and so on.
%
%       plt_flag : a flag to indicate whether or not to plot the solution.
%               == 0 - no plots
%               != 0 - turn on the plots
%
%       Nrefine  : number of times the solution was refined.
%
% See also: get_u2d, construct_pre_post_process_errors2d and postprocess2d.

clear all;

% ---------------------------------------------------------------------------- %
% Input paramters
% ---------------------------------------------------------------------------- %

%runname = 'SL4-Short';   % Used in name of output file (see fname below)
runname  = '_';           % Used in name of output file (see fname below)
plt_flag = 0;             % plotting flag
Nrefine  = 7;             % Number of times to refine

% ---------------------------------------------------------------------------- %
% Pull the order and starting resolution, and computational domain for this convergence study.
% ---------------------------------------------------------------------------- %
IniFileName = ['output', runname, num2str(0, '%03d'), '/parameters.ini'];
INI         = ConvertIniFile2Struct(IniFileName);                % Struct of structs
mpp         = ( str2num(INI.dogparams.space_order) )-1;     % Polynomial order
meqn        = ( str2num(INI.dogparams.meqn) );              % Number of equations
kmax        = (mpp+1)^2;

% computational domain
xbounds     = [str2num(INI.grid.xlow), str2num(INI.grid.xhigh)];
ybounds     = [str2num(INI.grid.ylow), str2num(INI.grid.yhigh)];
ElmtsStart  = [str2num(INI.grid.mx),   str2num(INI.grid.my   )]; % Starting grid resolution

% ---------------------------------------------------------------------------- %
% Main loop
% ---------------------------------------------------------------------------- %
for M=1:Nrefine

    Elmts = [ElmtsStart(1), ElmtsStart(2)]*2^(M-1);     % Number of grid points
    mx = Elmts(1); my = Elmts(2);

    Nold  = [ElmtsStart(1), ElmtsStart(2)]*2^(M-2);     % old number of grid points (for convergence study)

    % Vector of unkowns (Read in the DoGPack data file)
    fname = ['output', runname, num2str(M-1, '%03d') '/q0001.dat'];

    % Read the state variables
    [u0, time0]  = get_u2d( ['output', runname, num2str(M-1, '%03d')], 0, Elmts, meqn, mpp+1 );
    [u,  time ]  = get_u2d( ['output', runname, num2str(M-1, '%03d')], 1, Elmts, meqn, mpp+1 );

    % Evaluate the basis functions, and compute errors for the pre- and post-processed solution
    [XX,YY,   Approx,   Error,   er2,  erinf  ] = construct_pre_post_process_errors2d( ...
        u, time, mpp, xbounds, ybounds );
    [XX,YY, ppApprox, ppError, er2pp,  erinfpp, uex] = postprocess2d( ...
        u, time, mpp, xbounds, ybounds   );

    projection_error_l2 = norm( reshape(u-u0,mx*my*kmax*meqn,1), 2 ) / norm( reshape(u0, mx*my*kmax*meqn, 1 ), 2 );

    % Extract convergence rate based on the iteration number
    if( M==1 )
        order_projection_l2 = 0;
        pporderl2   = 0;        orderl2   = 0;
        pporderlinf = 0;        orderlinf = 0;
    else
        order_projection_l2 = log( old_projection_error / projection_error_l2)/log(Elmts(1)/Nold(1));
        pporderl2   = log(oldl2pp/er2pp        )/log(Elmts(1)/Nold(1));
        pporderlinf = log(oldlinfpp/erinfpp    )/log(Elmts(1)/Nold(1));
        orderl2     = log(oldl2/er2            )/log(Elmts(1)/Nold(1));
        orderlinf   = log(oldlinf/erinf        )/log(Elmts(1)/Nold(1));
    end

    fprintf('* -- (mx, my) = (%d,%d) -- *\n', Elmts(1), Elmts(2) )
    fprintf('Post-processed L-2 error         = %2.10e;     order=%2.3f\n', er2pp,   pporderl2  )
    fprintf('Post-processed L-inf error       = %2.10e;     order=%2.3f\n', erinfpp, pporderlinf)
    fprintf('\n')

    fprintf('Pre-Post-processed L-2 error     = %2.10e;     order=%2.3f\n', er2,   orderl2  )
    fprintf('Pre-Post-processed L-inf error   = %2.10e;     order=%2.3f\n', erinf, orderlinf)
    fprintf('L-2 error (with I.C. projection) = %2.10e;     order=%2.3f\n', projection_error_l2, order_projection_l2 )
    fprintf('* ---------------- *\n\n')

    % Save the old errors (for convergence study)
    oldl2pp   = er2pp;    oldl2   = er2;
    oldlinfpp = erinfpp;  oldlinf = erinf;
    old_projection_error = projection_error_l2;

    % ---------------------------------------------------------------------------- %
    % Plot stuff, if requested
    % ---------------------------------------------------------------------------- %
    if( plt_flag )

        % Plot the solution, pre-postprocessing
        figure(1);  clf;
        plot( Y, Approx, 'go', 'LineWidth', 3 );
        hold on;
        plot( Y, q_init1d( Y-time ), 'r-', 'LineWidth', 3 );
        hold off;

        % Plot the post-processed errors
        figure(2)
        hold on;
        if M==1
            clf;
            semilogy(Ypp, ppError)
        elseif M==2
            semilogy(Ypp, ppError, 'r')
        elseif M==3
            semilogy(Ypp, ppError, 'g')
        elseif M==4
            semilogy(Ypp, ppError, 'c')
        elseif M==5
            semilogy(Ypp, ppError, 'b')
        end
        hold off;
        xlabel('X')
        ylabel('|Post-Processed Error|')
        title('Pointwise Post-Processed Error')
        axis( [xbounds(1) xbounds(2) 1e-8 1] );

        % Plot the pre-post-processed errors
        figure(3)
        hold on;
        if M==1
            clf;
            semilogy(Y, Error)
        elseif M==2
            semilogy(Y, Error, 'r')
        elseif M==3
            semilogy(Y, Error, 'g')
        elseif M==4
            semilogy(Y, Error, 'c')
        elseif M==5
            semilogy(Y, Error, 'b')
            legend(['N=',num2str(ElmtsStart,'%d')], ...
                ['N=',num2str(2*ElmtsStart,'%d')], ...
                ['N=',num2str(4*ElmtsStart,'%d')], ...
                ['N=',num2str(8*ElmtsStart,'%d')], ...
                ['N=',num2str(16*ElmtsStart,'%d')] );
        end
        hold off;
        xlabel('X')
        ylabel('|Pre-Post-processed Error|')
        title('Pointwise Pre-Post-Processed Error')
        axis( [xbounds(1) xbounds(2) 1e-8 1] );

        % Plot the modification errors
        Diff   = Approx-ppApprox;
        figure(4)
        hold on;
        if M==1
            clf;
            plot(Ypp, Diff)
        elseif M==2
            plot(Ypp, Diff, 'r')
        elseif M==3
            plot(Ypp, Diff, 'g')
        elseif M==4
            plot(Ypp, Diff, 'c')
        elseif M==5
            plot(Ypp, Diff, 'b')
        %   legend('N=8', 'N=16', 'N=32','N=64','N=128');
            legend(['N=',num2str(ElmtsStart,'%d')], ...
                ['N=',num2str(2*ElmtsStart,'%d')], ...
                ['N=',num2str(4*ElmtsStart,'%d')], ...
                ['N=',num2str(8*ElmtsStart,'%d')], ...
                ['N=',num2str(16*ElmtsStart,'%d')] );
        end
        hold off;
        xlabel('X')
        ylabel('|u_h-u^*|')
        title('Modification errors')

    end % End of plotting loop

end % End of Nrefine loop

if( plt_flag )
    % Save the pretty pictures!
    print(2, '-depsc', [runname, '_Pointwise post-processed-errors.eps'] );
    print(3, '-depsc', [runname, '_Pointwise  pre-processed-errors.eps'] );
    print(4, '-depsc', [runname, '_Modification errors.eps']             );
end


figure(1);
clf
surf(XX, YY, ppApprox(:,:,1) );
%colormap('jet');
axis on; box on; grid off;
axis('equal');
axis([-0.04 1.04 -0.04 1.04]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['Post-process solution at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 1]);
shading flat;

figure(2);
clf
surf(XX, YY, ppError(:,:,1) );
colormap('jet');
axis on; box on; grid off;
axis('equal');
%axis([-0.04 1.04 -0.04 1.04]);
axis auto;
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['Errors in post-processed solution at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 1]);
shading flat;

figure(3);
clf
surf(XX, YY, Approx - ppApprox );
colormap('jet');
axis on; box on; grid off;
axis('equal');
%axis([-0.04 1.04 -0.04 1.04]);
axis auto;
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.25:2);
set(gca,'fontsize',16);
t1 = title(['Modification errors at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
caxis([0 1]);
shading flat;

