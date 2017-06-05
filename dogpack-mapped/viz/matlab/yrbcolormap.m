%
% yellow-red-blue color map 
%
c01 = 0:.1:1;
c10 = 1:-.1:0;
c0 = 0*c01;
c1 = 1 + c0;
yrb = [c1' c10' c0'; c10' c0' c01'];
colormap(yrb)

% As of October 3, 2013, applications that call this file include:
%./apps/2d/euler/mach_reflection/plotq2_cart.m:yrbcolormap
%./apps/2d/euler/radial_shock/plotq2_cart.m:yrbcolormap
%./apps/2d/euler/riemann_test/plotq2_cart.m:yrbcolormap
%./apps/2d/euler/test_exact/plotq2_cart.m:  yrbcolormap
%./apps/2d/shallow_water/riemann_test/plotq2_cart.m:yrbcolormap
%./apps/plasma/2d/twofluid/g05/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/g05/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/g05/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/g05/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/g05/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/g05/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/g10/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/i10e5/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/p05/plotq2_cart.m:yrbcolormap;
%./apps/plasma/2d/twofluid/p10/plotq2_cart.m:yrbcolormap;
%./viz/matlab/plotq2_cart.m:yrbcolormap
%
% I think it would be nice to move this outside of the top-level matlab
% library.  For example, we create a 'color options', or perhaps a 'plot
% options' directory, because this doesn't really have anything to do with the
% top-level functions, nor the size and shape of the polynomials.  (-DS)
