function graphrecon(problem)

if(~exist('problem','var'))
  problem='rescaled'
end

switch problem
  case 'GEM';
    timemax=50;
    ratemax=.16;
  case 'rescaled';
    timemax=80;
    ratemax=.15;
  otherwise
    error('no such problem');
end

switch problem
  case 'GEM';
    a=[[32  , .14 , 30]; ...
       [32.5, .13 , 10]; ...
       [33  , .13 , 6 ]; ...
       [33.5, .14 , 3 ]; ...
       [36  , .13 , 1 ]; ...
       [37.5, .14 , .3]]; ...
    xticklabel = {30, 10, 6, 3, 1, .3};
  otherwise
    a=[[ 68,  .06, 200000000000]; ...
       [ 64,  .07,  200]; ...
       [ 58,  .075,  60]; ...
       [ 53,  .08,   20]; ...
       [ 50,  .09,    6]; ...
       [ 51,  .08,    2]; ...
       [ 52,  .06,   .6]; ...
       [ 43,  .09,   .2]; ...
       [ 40,  .125, .06]; ...
       [ 38,  .13,  .02]; ...
       [ 35,  .13,  0.000000000001]];
    xticklabel = {'infinity     ','200  ','60','20','6','2','.6','.2','.06 ',' .02','0'};
  end
  
  rates = rescale(2./a(:,3),0.23);
  clf;
  figure(1);
  %plot(rates,1.662./a(:,1),'--rs','linewidth',2);
  plot(rates,a(:,1),'-bs','linewidth',3,'MarkerSize',20);
  ylim([0,timemax]);
  set(gca,'XTick',rates);
  set(gca,'XTickLabel',xticklabel);
  title('time until 30% of flux is reconnected',...
    'fontsize',26,'fontweight','bold');
  xlabel('isotropization period per gyroperiod',...
    'fontsize',20,'fontweight','bold');
  ylabel('time per gyroperiod',...
    'fontsize',20,'fontweight','bold');
  set(gca,'fontsize',20,'fontweight','bold','linewidth',3);
  set(gca,'YTick',0:10:80);
  limits = axis();
  xpos = limits(1)+.02*(limits(2)-limits(1));
  ypos = limits(3)+.09*(limits(4)-limits(3));
  text(xpos, ypos, ['(' problem ' problem)'],'fontsize',20,'fontweight','bold');
  print('-depsc', [getenv('HOME') '/figures/recon_time']);

  figure(2);
  plot(rates,a(:,2),'-bs','linewidth',3,'MarkerSize',20);
  ylim([0,ratemax]);
  set(gca,'XTick',rates);
  set(gca,'XTickLabel',xticklabel);
  title('peak rate of reconnection',...
    'fontsize',26,'fontweight','bold');
  xlabel('isotropization period per gyroperiod',...
    'fontsize',20,'fontweight','bold');
  ylabel('peak electric field (per B_0v_A)',...
    'fontsize',20,'fontweight','bold');
  set(gca,'fontsize',20,'fontweight','bold','linewidth',3);
  set(gca,'YTick',0:.02:.16);
  limits = axis();
  xpos = limits(1)+.02*(limits(2)-limits(1));
  ypos = limits(3)+.09*(limits(4)-limits(3));
  text(xpos, ypos, ['(' problem ' problem)'],'fontsize',20,'fontweight','bold');
  print('-depsc', [getenv('HOME') '/figures/peak_recon_rate']);

end

function out=rescale(a,b)
  % tanh(b*ln(a))
  ab = a.^b;
  ab2 = ab.^2;
  out=(ab2-1)./(ab2+1);
end
