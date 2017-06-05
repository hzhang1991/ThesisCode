%
%  colormap (cold -> hot)
%
z16 = zeros(1,16);  o32 = 32*ones(1,32);
hj  = [[z16,(0:2:32),o32];[(0:1:32),(31:-1:0)];[o32,(32:-2:0),z16]]'/32;

colormap(hj);