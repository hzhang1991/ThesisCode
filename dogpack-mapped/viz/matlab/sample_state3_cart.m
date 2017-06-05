function qvals = sample_state3_cart(q, method_order,kmax,LegVals)
%SAMPLE_STATE3_CART.  Convert coefficients to point values.
%
% 

   mx = size(q,1);
   my = size(q,2);
   mz = size(q,3);
   meqn = size(q,5);  
   mpts = size(LegVals,2);   
   points_per_dir = round((mpts)^(1/3));   
   qvals = zeros(mx*my*mz*mpts,meqn);   
   index = 0;  
   
   for k=1:mz
       for m1=1:points_per_dir
           for j=1:my
               for m2=1:points_per_dir
                   for i=1:mx
                       for m3=1:points_per_dir
                           
                           m = m3 + points_per_dir*(m2-1 + points_per_dir*(m1-1));
                           
                           index = index+1;
                           for n=1:meqn
             
                               qvals(index,n) = 0;
                               
                               for ell=1:kmax
                                   qvals(index,n) = qvals(index,n) + q(i,j,k,ell,n)*LegVals(ell,m);
                               end  
                               
                           end
                       end
                   end
               end
           end
       end
   end
           
   qvals = reshape(qvals,mx*points_per_dir,my*points_per_dir,mz*points_per_dir,n);

end
