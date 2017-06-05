function qvals = sample_state4_hybrid(q, method_order, kmax, LegVals )
%SAMPLE_STATE4_HYBRID   Sample the basis functions.
%
% LegVals need to be constructed for each triangle sampled.

   mx     = size(q,1);
   my     = size(q,2);
   melems = size(q,3);
   meqn   = size(q,4);

   mpts   = size(LegVals,2);
   qvals = zeros(melems*mpts, meqn);

   index = 0;
   for m1=1:mx
   for m2=1:my
   for i=1:melems
     for mp=1:mpts
       index = index+1;
       for me=1:meqn
         qvals(index,me) = 0;
         for k=1:kmax
           qvals(index,me) = qvals(index,me) + q(m1,m2,i,me,k)*LegVals(k,mp);
         end  
       end
     end
   end
   end
   end

end
