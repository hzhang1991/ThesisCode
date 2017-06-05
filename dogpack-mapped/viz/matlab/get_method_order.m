function method_order = get_method_order(kmax, ndims)
%GET_METHOD_ORDER.   Get the method order.
%
% method_order = GET_METHOD_ORDER(kmax, ndims) returns the order of the method
% by counting the number of polynomials kmax used in the Galerkin expansion.  
% This is the inverse function of GET_KMAX.
%
% Input:
%
%     kmax  - number of polynomials used in an expansion.
%     ndims - number of dimensions the expansion was written in.  Currently
%             only supports ndims = 1,2 or 3.
%
% Output:
%
%     method_order - the order of the method.
%
% See also: GET_KMAX.

    if( ndims==1 )
        method_order = kmax;
    elseif( ndims==2 )
        method_order = round((-1+sqrt(1+8*kmax))/2);
    elseif( ndims==3 )
        dd = (81*kmax + 3*sqrt(729*kmax^2-3))^(1/3);
        method_order = round(dd/3 + 1/dd - 1);
    else
        error(['unsupported ndims:' num2str(ndims)]);        
    end
  
end
