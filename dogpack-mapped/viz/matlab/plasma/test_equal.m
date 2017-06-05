% This belongs in the main matlab directory
%
% test if two floats are equal
function bool=test_equal(a,b)
  eps=1e-12;
  mag = abs(a)+abs(b);
  diff = abs(a-b);
  if(mag<1)
    if(diff<eps)
      bool=true;
    else
      bool=false;
    end
  else
    if(diff/mag<eps)
      bool=true;
    else
      bool=false;
    end
  end
end
