function ke = fit_coral_iCA(time, Obsdata, p)
global FIG_NUM;

kvar_init(1) = 10;         % /s, forward CO2 hydration rate constant catalyzed by C

%scale parameters
kvar_init = kvar_init .* p.scale;

mink(1) = kvar_init(1)./1001;
maxk(1) = kvar_init(1).*1001;

%set options for lsqcurvefit
      opts = optimset('lsqcurvefit');      
      opts = optimset(opts,'DiffMaxChange',100);
    % opts = optimset(opts,'DiffMinChange',5E-7);
      opts = optimset(opts,'Display','iter');
      opts = optimset(opts,'MaxIter',200); 
      opts = optimset(opts,'Diagnostics','off');
      opts = optimset(opts,'MaxFunEvals',8000);   
      opts = optimset(opts,'TolFun',1e-22);
      opts = optimset(opts,'TolX',1E-10);
      
[kfit, resnorm, residual, exitflag, output, lambda, Jac] = lsqcurvefit('coral_iCA_lsq', kvar_init, time, Obsdata, mink, maxk, opts, p);

[m,n] = size(residual);
residual = reshape(residual,m*n,1);
FIG_NUM = FIG_NUM + 1;
figure(FIG_NUM)
hist(residual,20),title('Residuals');

%estimate 95% confidence intervals on fitted parameters
ci = nlparci(kfit, residual,'jacobian', Jac);
sd = kfit' - ci(:,1);
kfit = kfit./p.scale;
sd = sd./p.scale';

ke = cat(2,kfit',sd);
return