% meqn: number of variables
% mx: number of mesh points in x dimension
% my: number of mesh points in y dimension
% method1: order of accuracy in space
%
function [qfull,time]=readq2(n1,outputdir_in,datafmt,...
  meqn,mx,my,method1,components)

  %if(nargin<9)
  %  READ_FROM_RESTART=1;
  %end

  READ_FROM_RESTART=1;

  %basefilename = [outputdir sprintf('/q%04d',n1)];
  %[q,time]=read_state2_cart(datafmt, basefilename, '/q',...
  %  mx, my, meqn, get_kmax(method1), components);
  [q,time]=read_state2_cart(datafmt, outputdir, n1, 'q',...
    mx, my, meqn, get_kmax(method1), components);
  qfull = sample_state2(q, method1);

  global outputdir
  if(nargin>=2)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end
  if(nargin<7)

     error('This code is not tested.');

     %fids = fopen([outputdir '/qhelp.dat'],'r');
     %file_meqn = fscanf(fids,'%d',1);
     %file_nplot = fscanf(fids,'%d',1);
     %file_meth1 = fscanf(fids,'%d',1);
     %file_mx    = fscanf(fids,'%d',1);
     %file_my    = fscanf(fids,'%d',1);
     %file_datafmt = 1;
     %file_datafmt = fscanf(fids,'%d',1);
     %fclose(fids);
     [file_mx, ...
      file_my, ...
      file_xlow, ...
      file_xhigh, ...
      file_ylow, ...
      file_yhigh, ...
      file_method1, ...
      file_meqn, ...
      file_nout, ...
      file_mx_out, ...
      file_my_out, ...
      file_datafmt]...
          =get_plot_params2(outputdir);

     if(nargin<7); method1 = file_method(1); end;
     if(nargin<6); my      = file_my; end;
     if(nargin<5); mx      = file_mx; end;
     if(nargin<4); meqn    = file_meqn; end;
     if(nargin<3); datafmt = file_datafmt; end;
  end
  if(nargin<8)
    components=1:meqn;
  end
  basefilename = [outputdir sprintf('/q%04d',n1)];

  if(READ_FROM_RESTART)
    tic;
    [q,time]=read_state2(datafmt, basefilename, '/q',
      mx, my, meqn, method1, components);
    toc;
    tic;
    qfull = sample_state2(q, method1);
    toc;
  else 
  % 
    tic;
    [qfull,time]=read_output(datafmt, basefilename, ...
      mx, my, meqn, method1, components);
    toc;
  end
  %size(qfull)
end
    
