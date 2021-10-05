function CO2_out = uncat_lsq(kuf, time, par)
% computes Cspecies time course for given kuf
% using analytical solution of Silverman 1974

time = time - time(1,1);                                %reset time to start at zero
CO2tot = [1 1 1] * par.cinit;                           %total CO2ext conc at start of phase
Tau = ([0 0.5 1] * par.cinit)./ CO2tot;                 %O18 atom fraction in CO2 
c = ([0 0 1] * par.cinit)./CO2tot;                      %fraction of doubly labelled CO2

Hp = 10^-par.pHe;                                       % H+ concentration
kur = kuf .* Hp./par.K1;                                %reverse rate const, HCO3-> CO2           
f = 1/((Hp/par.K1) +  1 + (par.K2/Hp));                 %fraction of dic as HCO3- (truely, not HCO3- + CO32-)
ktau = -kur .* f .* (1/3);                              %decay constant for Tau
kc = -kur .* f .* (2/3);                                %decay constant for c

Tau_t = Tau .* exp(ktau.*time);                         %decay of tau over time
c_t = c .* exp(kc.*time);                               %decay of c over time

CO2_out(3,:)= c_t .* CO2tot;                            %CO2 mass 49
CO2_out(2,:) = 2.*(Tau_t.*CO2tot - CO2_out(3,:));       %CO2 mass 47
CO2_out(1,:) = CO2tot - CO2_out(2,:) - CO2_out(3,:);    %CO2 mass 45    

return