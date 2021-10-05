function uncat_vals = fit_uncat(time, Obsdata, par)
global FIG_NUM;
%fits 18O isotope exchange data to CO2 +H20 -> HCO3 + H equation
%and derives calibration (or response) factors for the CO2 signals

k0(1) = 5.72E-2;                          %initial guess for rate constant for HCO3 formation from CO2 + H20, 
%k0(2) = (10^-pH1)*2.66E4 ;               %/s ;rate constant for CO2 formation from HCO3 by reaction w/ H+, from Zeebe and Wolf-Gladrow

mink(1) = k0(1)/10000;
maxk(1) = k0(1)*10000;

%set options for lsqcurvefit
      opts = optimset('lsqcurvefit');      
      opts = optimset(opts,'Display','iter');
      opts = optimset(opts,'MaxIter',500); 
      opts = optimset(opts,'Diagnostics','off');
      opts = optimset(opts,'MaxFunEvals',24000);   
      opts = optimset(opts,'TolFun',1e-30);
%      opts = optimset(opts,'TolX',[0.001 0.001 0.1 0.01]);
      
[kfit, resnorm, residual, exitflag, output] = lsqcurvefit('uncat_lsq', k0 ,time, Obsdata, mink, maxk, opts, par);
%fit data with single exchange rate
CO2pred = uncat_lsq(kfit(1), time ,par);

C45pred = CO2pred(1,:);
C47pred = CO2pred(2,:);
C49pred = CO2pred(3,:);

FIG_NUM = FIG_NUM+1;
figure(FIG_NUM);
plot(time,Obsdata(1,:),'bo',time,C45pred, time, Obsdata(2,:), 'go', time,C47pred,'g', time, Obsdata(3,:), 'ro', time, C49pred,'r'), xlabel('time'), ylabel('CO2 (mol/cm3)');

uncat_vals = kfit;
return
