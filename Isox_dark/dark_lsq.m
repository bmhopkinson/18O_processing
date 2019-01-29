function CO2out = dark_lsq(kvar, time, par)
%interface from lsq to ode call; set up problem to call ode solver
global FIG_NUM;
%testing svn
x = 1;

kvar = kvar ./ par.scale;       %unscale parameters
par.kcf = kvar(1);              %/s, forward CO2 hydration rate constant catalyzed by CA inside the cell, initial guess
par.fc = kvar(2);               %/s, rate constant for diffusive CO2 flux
par.fb = kvar(3);               %/s, rate constant for diffusive HCO3 flux

CO2tot = [1 1 1] * par.cinit(1:3,1);           %total CO2ext conc at start of phase
Cinit(1:3,1) = par.cinit;
%calculate HCO3 concentrations assuming chemical and isotopic (apply binomial) equilibrium among species, but not w/ current CO2;
Htot = par.DIC - CO2tot;                                        %HCO3- + CO32- concentration, call HCO3- and treat the same, making corrections to hyd/dehyd rates to account for it.
Cinit(4,1) = Htot * (1 - par.Taub)^3;                            %H61, row 4
Cinit(5,1) = Htot * (3 * (par.Taub) * (1-par.Taub)^2);            %H63, row 5
Cinit(6,1) = Htot * (3 * (par.Taub^2) * (1 - par.Taub));          %H65, row 6
Cinit(7,1) = Htot * (par.Taub^3);                                %H67, row 7

%approximate initial CO2in and HCO3in as equivalant to external
Cinit(8:10,1) = Cinit(1:3,1);
Cinit(11:14,1) = Cinit(4:7,1);

%use ode solver to deterimine time course of isotope exchange
options = odeset('RelTol', 1E-4, 'AbsTol', 1E-9);         %set ode options to obtain smooth (non-oscillating or diverging) solution
[t_ode, Cspecies] = ode15s(@darkderiv, time, Cinit, options, par);
Cspecies = Cspecies';                  %ode gives one row per time point, transpose to return to normal structure: each row is the time series of a single Cspecies

  
if par.return == 0
CO2out = Cspecies(1:3,:);

elseif par.return == 1     
CO2out = Cspecies;
%write predicted data out to file, this is the fitted value
outfile = 'iso_exch_dark.out';
fout_id = fopen(outfile,'w');
fprintf(fout_id,'time\t C45ext\t C47ext\t C49ext\t H61ext\t H63ext\t H65ext\t H67ext\t C45in\t C47in\t C49in\t H61in\t H63in\t H65in\t H67in\n');
reps = size(Cspecies,2);
for i = 1:reps
    fprintf(fout_id,'%6.2f\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\n',time(i),Cspecies(:,i));
end
fclose(fout_id);      

end

return
