function Cfit = kap_lsq(kcfit, time, par)

%ode solver requires initial conditions as a single vector
Cinit(1:15,1) = [par.ce_init; par.be_init; par.ci_init; par.bi_init; par.O2_init];

%tranfer fitting parameter to par structure
%fprintf(1,'kcfit %e\n',kcfit);
par.kcf_fit = kcfit;
fb =(1/ (1+ par.K2./(10^-par.pHe)));                    % fraction of "h" pool that is HCO3- to correct HCO3 dehydration rates. 
par.kcr_fit = par.kcf_fit * fb * (10^-par.pHe)./par.K1;              % catalyzed rate of HCO3 dehydration

options = odeset('RelTol', 1E-4, 'AbsTol', 1E-8);         %set ode options to obtain smooth (non-oscillating or diverging) solution

[t_ode, Cspecies] = ode15s(@kapderiv, time, Cinit, options, par);
Cspecies = Cspecies';                                                       %ode gives one row per time point, transpose to return to normal structure: each row is the time series of a single Cspecies


if par.flag == 0                  %format required for lsqcurvefit
Cfit = Cspecies(1:3,:);
elseif par.flag ==1               %full Cspecies data and ode time
Cfit(1,:)   = t_ode;
Cfit(2:16,:) = Cspecies;
end

return