function qvals = sample_state2_cart(q, method_order, kmax, LegVals)
%SAMPLE_STATE2_CART.  Convert the 2D coefficients into point values.
%
% The purpose of this call is to convert Galerkin coefficients into point
% values.  This function turns coefficients into point values on 2D Cartesian
% grid.  This is the primary routine called by PLOTDOG2.
%
% This routine assumes that the user wants to plot the same quadrature points 
% in each element in the entire domain.
%
% Input:
%
%     q - state variable.  size(q) = [mx, my, kmax, meqn].
%
%     method_order - order of the method.  (Currently not used in this
%                   function.)
%
%     kmax - number of polynomials in q.  (size(q,3) == kmax).
%
%     LegVals - Legendre polynomials evaluated at desired points.
%
% Output:
%
%     qvals - pointwise values of state vector.  
%             size( qvals ) = [ mx*points_per_dir, my*points_per_dir, meqn ]
%
%             See the for loops to see the order in which these points are
%             stacked.
%
% See also: PLOTDOG2, READ_STATE2_CART, SAMPLE_STATE2_CART and SAMPLE_STATE_CART2

   mx   = size(q,1);
   my   = size(q,2);
   meqn = size(q,4);  
   mpts = size(LegVals,2);   
   points_per_dir = sqrt(mpts);   
   qvals = zeros(mx*my*mpts,meqn);   
   index = 0;  
   
   for j=1:my
     for m1=1:points_per_dir
       for i=1:mx
         for m2=1:points_per_dir
           
           m = m2 + points_per_dir*(m1-1);
           index = index+1;
           for n=1:meqn
             
             qvals(index,n) = 0;
             for k=1:kmax
               qvals(index,n) = qvals(index,n) + q(i,j,k,n)*LegVals(k,m);
             end  
             
           end
         end
       end
     end
   end

   % Reshape q for cleaner indexing:
   qvals = reshape(qvals, mx*points_per_dir, my*points_per_dir, meqn);

end
