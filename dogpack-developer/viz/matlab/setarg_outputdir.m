% call this batch file to ensure that the outputdir
% is set to outputdir_in if provided, and to the
% default value of 'output' if there is nothing else
% to fall back on.

global outputdir;
if(exist('outputdir_in','var'))
outputdir=outputdir_in;
elseif(isempty(outputdir))
outputdir='output'
end
if(~isdir(outputdir))
error(['no such directory: ' outputdir]);
end

% As of October 3, 2013, current applications that use this file include:
%
% ./apps/plasma/2d/twofluid/g05/matlab/plotout.m:  setarg_outputdir;
% ./apps/plasma/2d/twofluid/g10/matlab/plotout.m:  setarg_outputdir;
% ./apps/plasma/2d/twofluid/i10e5/matlab/plotout.m:  setarg_outputdir;
% ./apps/plasma/2d/twofluid/p05/matlab/plotout.m:  setarg_outputdir;
% ./apps/plasma/2d/twofluid/p10/matlab/plotout.m:  setarg_outputdir;
% ./viz/matlab/plasma/@GEMplotter/GEMplotter.m:  setarg_outputdir;
% ./viz/matlab/plasma/@Xplotter/plot_recon.m:  setarg_outputdir;
% ./viz/matlab/plasma/@Xplotter/plot_xpoint.m:  setarg_outputdir;
% ./viz/matlab/plasma/make_outputdirname.m:  setarg_outputdir
% ./viz/matlab/plasma/plot_xpress.m:  setarg_outputdir
% ./viz/matlab/plasma/plotall.m:  setarg_outputdir;
% ./viz/matlab/plasma/xpoint_state_to_vals.m:  setarg_outputdir;
