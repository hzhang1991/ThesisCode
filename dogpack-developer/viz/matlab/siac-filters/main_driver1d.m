%MAIN_DRIVER.
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
% See also: get_u1d, construct_pre_post_process_errors1d and postprocess1d.

clear all;

% ---------------------------------------------------------------------------- %
% Input paramters
% ---------------------------------------------------------------------------- %

%runname = 'SL4-Short';   % Used in name of output file (see fname below)
runname  = '';            % Used in name of output file (see fname below)
plt_flag = 1;             % plotting flag
Nrefine  = 6;             % Number of times to refine

% ---------------------------------------------------------------------------- %
% Pull the order and starting resolution, and computational domain for this convergence study.
% ---------------------------------------------------------------------------- %
IniFileName = ['output', runname, num2str(0, '%03d'), '/parameters.ini'];
INI         = ConvertIniFile2Struct(IniFileName);       % Struct of structs
mpp         = str2num(INI.dogparams.space_order)-1;     % Polynomial order

% computational domain
xbounds     = [str2num(INI.grid.xlow), str2num(INI.grid.xhigh)];
ElmtsStart  = str2num(INI.grid.mx);                     % Starting grid resolution

% ---------------------------------------------------------------------------- %
% Main loop
% ---------------------------------------------------------------------------- %
for M=1:Nrefine

    Elmts = ElmtsStart*2^(M-1);     % Number of grid points
    Nold  = ElmtsStart*2^(M-2);     % old number of grid points (for convergence study)

    % Vector of unkowns (Read in the DoGPack data file)
    fname = ['output', runname, num2str(M-1, '%03d') '/q0001.dat'];
    [u0, time0]  = get_u1d( ['output', runname, num2str(M-1, '%03d') '/q0000.dat'], Elmts, mpp+1 );
    [u,  time ]  = get_u1d( ['output', runname, num2str(M-1, '%03d') '/q0001.dat'], Elmts, mpp+1 );

    % Evaluate the basis functions, and compute errors for the pre- and post-processed solution
    [Y,     Approx, Error,   er2,    erinf  ] = construct_pre_post_process_errors1d(u, time, mpp, xbounds, Elmts);
    [Ypp, ppApprox, ppError, er2pp,  erinfpp] = postprocess1d( u, time, mpp, xbounds, Elmts );

    projection_error_l2 = norm( u-u0, 2 ) / norm( u0, 2 );

    % Extract convergence rate based on the iteration number
    if( M==1 )
        order_projection_l2 = 0;
        pporderl2   = 0;        orderl2   = 0;
        pporderlinf = 0;        orderlinf = 0;
    else
        order_projection_l2 = log( old_projection_error / projection_error_l2)/log(Elmts/Nold);
        pporderl2   = log(oldl2pp/er2pp)/log(Elmts/Nold);
        pporderlinf = log(oldlinfpp/erinfpp)/log(Elmts/Nold);
        orderl2     = log(oldl2/er2)/log(Elmts/Nold);
        orderlinf   = log(oldlinf/erinf)/log(Elmts/Nold);
    end

    fprintf('* -- Elmts = %d -- *\n', Elmts);
    fprintf('Post-processed L-2 error         = %2.10e;     order=%2.3f\n', er2pp,   pporderl2  )
    fprintf('Post-processed L-inf error       = %2.10e;     order=%2.3f\n', erinfpp, pporderlinf)
    fprintf('\n')

    fprintf('Pre-Post-processed L-2 error     = %2.10e;     order=%2.3f\n', er2,   orderl2  )
    fprintf('Pre-Post-processed L-inf error   = %2.10e;     order=%2.3f\n', erinf, orderlinf)
    fprintf('L-2 error (with I.C. projection) = %2.10e;     order=%2.3f\n', projection_error_l2, order_projection_l2 )
    fprintf('* ---------------- *\n\n');

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

% Save the pretty pictures!
print(2, '-depsc', [runname, '_Pointwise post-processed-errors.eps'] );
print(3, '-depsc', [runname, '_Pointwise  pre-processed-errors.eps'] );
print(4, '-depsc', [runname, '_Modification errors.eps']             );

