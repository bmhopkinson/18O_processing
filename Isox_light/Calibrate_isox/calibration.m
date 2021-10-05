function CFs = calibration(CO2_signals, CO2_labtot, enrich_factor, uncat_samps)
global FIG_NUM

cf_init(1) = 1.0;                              % C45 response factor (RF), initial guess
%cf_init(2) = 1.0;                              % RF C47
%cf_init(3) = 1.0;                              % RF C49

mink(1) = cf_init(1)/1E6;
%mink(2) = cf_init(2)/10;
%mink(3) = cf_init(3)/10;

maxk(1) = cf_init(1)*1E6;
%maxk(2) = cf_init(2)*10;
%maxk(3) = cf_init(3)*10;

samps = size(CO2_signals, 2);

x(1,:) = 1:1:samps;
Exp_value(1,(1:uncat_samps)) = enrich_factor .* CO2_labtot;
Exp_value(1,(uncat_samps:samps)) = CO2_labtot;

params = CO2_signals;

%this should really be a linear multiple regression I think, revise if
%needed.
%set options for lsqcurvefit
      opts = optimset('lsqcurvefit');      
      opts = optimset(opts,'Display','iter');
      opts = optimset(opts,'MaxIter',4000); 
      opts = optimset(opts,'Diagnostics','off');
      opts = optimset(opts,'MaxFunEvals',24000);   
      opts = optimset(opts,'TolFun',1e-30);
      opts = optimset(opts,'TolX',[1E-6]);
      
[cf_fit, resnorm, residual, exitflag, output] = lsqcurvefit('apply_cfs', cf_init, x, Exp_value, mink, maxk, opts, params);

fprintf(1,'CFs %6.4E\n',cf_fit);

CO2_sum = apply_cfs(cf_fit, x, params);

FIG_NUM = FIG_NUM + 1;
figure(FIG_NUM);
plot(x,Exp_value, x, CO2_sum, 'g');

CFs = cf_fit;

return