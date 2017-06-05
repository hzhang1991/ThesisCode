function kmax = get_kmax(space_order, ndims)
%GET_KMAX    Number of DG polynomials.
%
% kmax = GET_KMAX( space_order, ndims ) returns the number of 
%        basis functions per element that is used DoGPack.
%
% Our implementation uses the minimal number of polynomials to achieve
% space_order order of accuracy.
%
% Input:
%
%   space_order - spatial order of accuracy
%   ndims       - number of dimensions
%
% Output:
%
%   kmax - number of basis functions used per element
%
% See also: plotdog1, plotdog2, plotdog3, plotdog4

  if(ndims==1)
    kmax=space_order;
  elseif(ndims==2)
    kmax=(space_order*(space_order+1))/2;
  elseif(ndims==3)
    kmax=(space_order*(space_order+1)*(space_order+2))/6;
  else
    fprintf(1,'Warning: are you sure you want ndims = %d?\n', ndims);
    deg  = space_order-1;
    kmax = nchoosek( ndims + deg, deg );
%   error(['unsupported ndims:' num2str(ndims)]);        
  end

end
