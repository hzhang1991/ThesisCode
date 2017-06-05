function plotdog1(points_per_dir, outputdir_in, point_type, qhelpname_in, qname_in, plotq1name_in)
%PLOTDOG1    Generic 1D plotting routine for DoGPack.
%
% PLOTDOG1( points_per_dir, outputdir_in, point_type  )
%
% points_per_dir = number of points used per cell for plotting.  Default = 1
%                  without any arguments supplied 
%
%   outputdir_in = location for output directory.  default = 'output'
%
%     point_type = 1:   uniform points on each element  ( Default )
%                = 2:   Gauss-Legendre points on each element
%
%   qhelpname_in = name of qhelp file used to get basic parameters
%
%       qname_in = name of variable filename (i.e., [qname]0000.dat, [qname]0001.dat, ...)
%
%  plotq1name_in = name of file for additional ploting options (i.e., filename is plotq1name.m)
%
% In almost every application in DoGPack, there is a single script, plotq1.m 
% that this function calls once per frame.  There is a default script located
% here in the library.  That script is where a
% user is expected to modify extra desired plotting stuff, for example setting
% axes, linecolors, exact solutions etc.
%
% See also: plotq1, plotdog2, plotdog3, plotdog4 for the other dimensional plotters.

    % Parse the input parameters (set default values for points_per_dir,
    % point_type and output folder name:
    parse_input_args;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  QHELP.DAT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fids  = fopen([outputdir,'/',qhelpname],'r');
    if fids==-1
        error(['File  ',outputdir,'/',qhelpname,'  not found.']);
    end
    ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    if (ndims~=1)
        error(['Incorrect dimension, ndims must be 1. ndims = ',num2str(ndims)]);
    end
    GridType = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    GridTmp = GridType(1:9);
    if (GridTmp=='Cartesian')
        GridType='Cartesian   ';
    else
        error(['Incorrect GridType, in 1D it must be Cartesian. GridType = ',GridType]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp(' ');
    disp(['        GridType = ',GridType]);
    disp(['  points_per_dir = ',num2str(points_per_dir)]);
    disp(['      point_type = ',num2str(point_type)]);
    disp(['       outputdir = ',outputdir]);
    disp(['       qhelpname = ',qhelpname]);
    disp(['           qname = ',qname]);
    disp(['      plotq1name = ',plotq1name]);
    
    disp(' ');

    % more information from qhelp.dat.  This is particular to the 1D solver.
    meqn    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    maux    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    nplot   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    meth1   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    datafmt = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    %
    mx      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    fclose(fids);

    % Grid information
    mx_old = mx;
    mx     = mx*points_per_dir;
    dx_old = (xhigh-xlow)/mx_old;
    dx     = (xhigh-xlow)/mx;
    
    % TODO - this needs to be modified ...
    if( point_type == 1 )

        xc     = transpose(linspace(xlow + 0.5*dx,     xhigh - 0.5*dx,     mx    ));
        xc_old = transpose(linspace(xlow + 0.5*dx_old, xhigh - 0.5*dx_old, mx_old));
        
        % Linearly spaced points
        s1d    = -1.0 + (2.0*(1:points_per_dir)-1.0)/points_per_dir;

    else
        
        sq3 = sqrt(3);
        sq5 = sqrt(5);
        sq7 = sqrt(7);
   
        % quadrautre points (TODO - add in 6th order story ... )
        if (points_per_dir==1)
            s1d = 0.0;
        elseif (points_per_dir==2)
            s1d = [-1.0/sq3; 1.0/sq3];
        elseif (points_per_dir==3)
            s1d = [-sq3/sq5; 0.0e0; sq3/sq5];
        elseif (points_per_dir==4)
            s1d = [-sqrt(3.0+sqrt(4.8))/sq7; -sqrt(3.0-sqrt(4.8))/sq7; ...
                   sqrt(3.0-sqrt(4.8))/sq7;  sqrt(3.0+sqrt(4.8))/sq7];
        elseif (points_per_dir==5)
            s1d = [-sqrt(5.0 + sqrt(40.0/7.0))/3.0; ...
                   -sqrt(5.0 - sqrt(40.0/7.0))/3.0; ...
                   0.0; ...
                   sqrt(5.0 - sqrt(40.0/7.0))/3.0; ...
                   sqrt(5.0 + sqrt(40.0/7.0))/3.0];
        end

        xc = zeros(mx_old*points_per_dir,1);
        
        xline = (xlow+dx_old/2):dx_old:(xhigh-dx_old/2);

        kk=1;
        for i=1:mx_old                 
            xc( kk:(kk+points_per_dir-1) ) = xline(i)+(dx_old/2)*s1d(:);
            kk=kk+points_per_dir;
        end

    end


    % Sample basis functions on mesh
    phi = GetCart1Legendre(meth1, s1d );

    q=-1;
    
    m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
                ' ) ? ']);
    disp(' ')
    if isempty(m)
        m=1;
    end

    kn = 0;

    n = 0;
    nf = 0;
    n1 = -1;

    while(nf~=-1)
        nf  = input([ ' Plot which frame ( 0 - ',num2str(nplot),...
                      ' ) [type -1 or q to quit] ? ']);
        if isempty(nf)
            n1 = n1 + 1;
            nf = 0;
        else
            n1 = nf;
        end
        if n1> nplot
            disp(' ');
            disp(' End of plots ');
            disp(' ');
            n1 = nplot;
        end
  
        if (nf~=-1)
            
            %% Solution -- q
            % solution should be found in file
            %     outputdir/q[n1].dat
            zz_qname = '';
            for zzq=1:length(qname)-1
                zz_qname(zzq)=qname(zzq);
            end
            fname = [outputdir,'/',zz_qname,num2str(n1+10000),'.dat'];
            
            % replace the 1000's digit by the letter q
            fname(length(outputdir)+length(zz_qname)+2) = qname(length(qname));
  
            fids = fopen(fname,'r');
            if fids==-1
                error(['File  ',fname,'  not found.']);
            end
            
            time = fscanf(fids,'%e',1);
            qtmp = fscanf(fids,'%e',[1,inf]);
            fclose(fids);
            qtmp = transpose(qtmp);
            qcoeffs  = reshape(qtmp, mx_old, meqn, meth1);
            clear qtmp;
            qsoln = zeros(mx,meqn);
            for i=1:mx_old
                for me=1:meqn
                    for ii=1:points_per_dir
                        v1(1:meth1,1) = phi(ii,:);
                        v2(1:meth1,1) = qcoeffs(i,me,:);
                        qsoln((i-1)*points_per_dir+ii,me) = transpose(v1)*v2;
                    end
                end
            end
            
            %% Aux variables -- aux
            if (maux>0)
                fname(length(outputdir)+2) = 'a';
                fids = fopen(fname,'r');
                if fids==-1
                    error(['File  ',fname,'  not found.']);
                end
                time = fscanf(fids,'%e',1);
                atmp = fscanf(fids,'%e',[1,inf]);
                fclose(fids);
                atmp = transpose(atmp);
                acoeffs  = reshape(atmp,mx_old,maux,meth1);
                clear atmp;
                aux = zeros(mx,maux);
                for i=1:mx_old
                    for ma=1:maux
                        for ii=1:points_per_dir              
                            v1(1:meth1,1) = phi(ii,:);
                            v2(1:meth1,1) = acoeffs(i,ma,:);
                            aux((i-1)*points_per_dir+ii,ma) = transpose(v1)*v2;
                        end
                    end
                end
            end
  
            % USER SUPPLIED FUNCTION
            plotq1;     
            %run(plotq1_varname);
            %run([plotq1name,'.m']);
        end
  
    end
    disp(' ')

end
