%% magnetic potential -- mp
fids  = fopen([outputdir,'/mhdhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/mhdhelp.dat  not found.']);
end
gamma  = fscanf(fids,'%e',1);
use_ct_method  = fscanf(fids,'%d',1);
use_ct_limiter = fscanf(fids,'%d',1);

if (use_ct_method==1)  
  mp_meth1 = meth1;%+1;
  mp_kmax = (mp_meth1)*(mp_meth1+1)/2;
  mp_LegVals = GetUnstLegendre(mp_kmax,points_per_dir,z);
  
  [mp_data,time] = read_state2_unst(datafmt, outputdir, n1, 'mp', ...
                                    NumElems, NumPhysElems, 1, mp_kmax, 1:1);
  mpsoln = sample_state2_unst(mp_data, mp_meth1, mp_kmax, mp_LegVals);
  clear mp_data;
end
