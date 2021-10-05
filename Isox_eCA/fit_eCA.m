function ke = fit_eCA(time, Obsdata, par)
global FIG_NUM;

%initial guess of parameter values
ksf_init = 	selfStart(time, Obsdata,par);      % use self starter function to find good initial guess for eCA activity

%minimum and maximum values for parameters
mink(1) = ksf_init(1)./1000;
maxk(1) = ksf_init(1).*1000;

%set options for lsqcurvefit, then call the optimization function
      opts = optimset('lsqcurvefit');      
      opts = optimset(opts,'DiffMaxChange',.1);
      opts = optimset(opts,'DiffMinChange',1E-7);
      opts = optimset(opts,'Display','iter');
      opts = optimset(opts,'MaxIter',100); 
      opts = optimset(opts,'Diagnostics','off');
      opts = optimset(opts,'MaxFunEvals',400);   
      opts = optimset(opts,'TolFun',1e-20);
      opts = optimset(opts,'TolX',1E-10);
      
[kfit, resnorm, residual, exitflag, output, lambda, Jac] = lsqcurvefit('eCA_lsq', ksf_init, time, Obsdata, mink, maxk, opts, par);

%estimate 95% confidence intervals on fitted parameters
ci = nlparci(kfit, residual,'jacobian', Jac);
ke = cat(2,kfit',ci);

return