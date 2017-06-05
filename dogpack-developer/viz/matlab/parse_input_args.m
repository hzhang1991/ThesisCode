%PARSE_INPUT_ARGS.    Single script for parsing input arguments.
%
% Each function plotdog1, plotdog2 and plotdog3 uses a common input format for
% its plotting routines.  The purpose of this function is to consolidate the
% input parameters passed into each function, which includes how many points
% per direction, points_per_dir, the location of the output directory,
% outputdir_in and the point type: point_type_in.
%
% Input parameters:
%
%   nargin - number of input arguments provided to plotdog[n].
%
%   points_per_dir - points per grid element.  Default = 1.
%
%   outputdir_in - string identifying the output directory.  This can be
%   relative or an absolute pathname.  Default = 'output'.
%
%   point_type    = 1:   uniform points on each element
%                 = 2:   Gauss-Legendre points on each element. Default: 1.
%
% Output:
%
%   Default values for each of these parameters, as stated above.
%
% See also: plotdog1, plotdog2, plotdog3, plotdog4

outputdir  = 'output';
qhelpname  = 'qhelp.dat';
qname      = 'q';
plotq1name = 'plotq1';

if(nargin<1)
    points_per_dir = 1;
end 
  
if (ischar(points_per_dir))
    points_per_dir = str2num(points_per_dir);
end

if(nargin>1)
    outputdir=outputdir_in;
elseif(isempty(outputdir))
    outputdir='output';
end

if points_per_dir<1
    points_per_dir = 1;
end

if (nargin>2)
    if(ischar(point_type))
        point_type = str2num(point_type);
    end
else
    point_type=1;
end

if point_type~=1 && point_type~=2
    point_type = 1;
end  
  
if point_type==2 && points_per_dir>5
    disp('too many points chosen, using 5 points instead');
    points_per_dir=5;
end
  
if(nargin>3)
    qhelpname=qhelpname_in;
elseif(isempty(qhelpname))
    qhelpname='qhelp.dat';
end
  
if(nargin>4)
    qname=qname_in;
elseif(isempty(qname))
    qname='q';
end
  
if(nargin>5)
    plotq1name=plotq1name_in;
elseif(isempty(plotq1name))
    plotq1name='plotq1';
end

format long e;
