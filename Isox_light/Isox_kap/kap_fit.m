function kcf_fit = kap_fit(time, Obsdata, par)

kcf_init(1) = 2*par.kif;         % /s, forward CO2 hydration rate constant catalyzed by CA, during photosynthesis
 
mink(1) = kcf_init(1)./100;
maxk(1) = kcf_init(1).*100;

%set options for lsqcurvefit
      opts = optimset('lsqcurvefit');      
      opts = optimset(opts,'DiffMaxChange',50);
      opts = optimset(opts,'DiffMinChange',25);
      opts = optimset(opts,'Display','iter');
      opts = optimset(opts,'MaxIter',400); 
      opts = optimset(opts,'Diagnostics','off');
      opts = optimset(opts,'MaxFunEvals',8000);   
      opts = optimset(opts,'TolFun',1e-26);
      opts = optimset(opts,'TolX',1E-12);
      
[kcf_fit, resnorm, residual, exitflag, output, lambda, Jac] = lsqcurvefit('kap_lsq', kcf_init, time, Obsdata, mink, maxk, opts, par);


%estimate 95% confidence intervals on fitted parameters
ci = nlparci(kcf_fit, residual,'jacobian', Jac);
sd = kcf_fit - ci(:,1);

kcf_fit = cat(2,kcf_fit,sd);

kcf_fit;
return

