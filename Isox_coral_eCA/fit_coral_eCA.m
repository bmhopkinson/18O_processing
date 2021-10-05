function ke = fit_coral_eCA(time, Obsdata, par)
global FIG_NUM;

kvar_init(1) = 1;             %cm/s eCA activity  
kvar_init(2) = 10;                %/s, forward CO2 hydration rate constant catalyzed by CA, had to /1000 to make of a similar order of magnitude as other parameters for fitting procedure (lsqcurvefit)
kvar_init(3) = 1E-2;              %cm/s, CO2 permeability  
kvar_init(4) = 1E-4;            %cm/s, HCO3- permeability

%scale parameters
kvar_init = kvar_init .* par.scale;

mink(1) = kvar_init(1)./100;
maxk(1) = kvar_init(1).*100;
mink(2) = kvar_init(2)./100;
maxk(2) = kvar_init(2).*100;
mink(3) = kvar_init(3)./100;
maxk(3) = kvar_init(3).*100;
mink(4) = kvar_init(4)./100;
maxk(4) = kvar_init(4).*100;

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
      
[kfit, resnorm, residual, exitflag, output, lambda, Jac] = lsqcurvefit('coral_eCA_lsq', kvar_init, time, Obsdata, mink, maxk, opts, par);

[m,n] = size(residual);
residual = reshape(residual,m*n,1);
FIG_NUM = FIG_NUM + 1;
figure(FIG_NUM)
hist(residual,20),title('Residuals');

%estimate 95% confidence intervals on fitted parameters
ci = nlparci(kfit, residual,'jacobian', Jac);
sd = kfit' - ci(:,1);
kfit = kfit./par.scale;
sd = sd./par.scale';

ke = cat(2,kfit',sd);
return