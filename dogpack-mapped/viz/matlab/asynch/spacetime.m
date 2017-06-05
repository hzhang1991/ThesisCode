function spacetime(outputdir_in)

  global outputdir
  outputdir='output';

  if(nargin==1)
    outputdir=outputdir_in;
  elseif(isempty(outputdir))
    outputdir='output';
  end

  format long e;

  disp(' ');
  disp(['       outputdir = ',outputdir]);
  disp(' ');
  
  fids  = fopen([outputdir,'/spacetime.dat'],'r');
  if fids==-1
    error(['File  ',outputdir,'/spacetime.dat  not found.']);
  end
  
  figure(1)
  clf;
  plot([0 1 1 0 0],[0 0 1 1 0],'k-');
  axis([0 1 0 1]);
  set(gca,'fontsize',16);
  t1=title('Spacetime evolution');
  set(t1,'fontsize',16);
  set(gca,'xtick',0:0.25:1);
  set(gca,'ytick',0:0.25:1);
  set(gca,'plotboxaspectratio',[1 1 1]);
  hold on;
  
  z = fscanf(fids,'%e',5);
  while ~feof(fids) 
    i  = z(1);
    t1 = z(2);
    t2 = z(3);
    x1 = z(4);
    x2 = z(5);
    if (i==41)
      p1=plot([x1 x2 x2 x1 x1],[t1 t1 t2 t2 t1],'r-'); 
    else
      p1=plot([x1 x2 x2 x1 x1],[t1 t1 t2 t2 t1],'b-');      
    end
    set(p1,'linewidth',2);
    pause(1.0e-12);
    z = fscanf(fids,'%e',5);
  end
    