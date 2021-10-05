function CO2out = eCA_lsq(kvar, time, par)
global FIG_NUM;
par.ksf = kvar(1);              % cm3/s, forward CO2 hydration rate constant catalyzed by eCA, initial guess

CO2tot = [1 1 1] * par.cinit(1:3,1);                            %total CO2ext conc at start of phase
Cinit(1:3,1) = par.cinit;
%calculate HCO3 concentrations assuming chemical and isotopic (apply binomial) equilibrium among species, but not w/ current CO2;
Htot = par.DIC - CO2tot;                                        %HCO3- + CO32- concentration, call HCO3- and treat the same, making corrections to hyd/dehyd rates to account for it.
Cinit(4,1) = Htot * (1 - par.Taub)^3;                           %B62, row 4
Cinit(5,1) = Htot * (3 * (par.Taub) * (1-par.Taub)^2);          %B64, row 5
Cinit(6,1) = Htot * (3 * (par.Taub^2) * (1 - par.Taub));        %B66, row 6
Cinit(7,1) = Htot * (par.Taub^3);                               %B68, row 7

%approximate initial CO2in and HCO3in as equivalant to external
Cinit(8:10,1)  = Cinit(1:3,1);       % CO2 surface (cs)
Cinit(11:14,1) = Cinit(4:7,1);       % HCO3- surface (bs)
Cinit(15:17,1) = Cinit(1:3,1);       % CO2 intracellular (ci)
Cinit(18:21,1) = Cinit(4:7,1);       % HCO3- intracellular (bi)

%use ode solver to deterimine time course of isotope exchange
options = odeset('RelTol', 1E-4, 'AbsTol', 1E-9);         %set ode options to obtain smooth (non-oscillating or diverging) solution
[t_ode, Cspecies] = ode15s(@eCAderiv, time, Cinit, options, par);
Cspecies = Cspecies';                  %ode gives one row per time point, transpose to return to normal structure: each row is the time series of a single Cspecies
  
if par.return == 0
  CO2out = Cspecies(1:3,:);

elseif par.return == 1     
  CO2out = Cspecies;

  %write predicted data out to file, this is the fitted value
  outfile = strcat(par.infile,'_eCA_fit.out');
  fout_id = fopen(outfile,'w');
  fprintf(fout_id,'time\t C45e\t C47e\t C49e\t B62e\t B64e\t B66e\t B68e\t C45s\t C47s\t C49s\t B62s\t B64s\t B66s\t B68s\t C45i\t C47i\t C49i\t B62i\t B64i\t B66i\t B68i\n');
  reps = size(Cspecies,2);
  for i = 1:reps
    fprintf(fout_id,'%6.2f\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\n',time(i),Cspecies(:,i));
  end
  fclose(fout_id);      
end

return
