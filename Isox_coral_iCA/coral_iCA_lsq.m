function CO2out = coral_iCA_lsq(kvar, time, p)
%interface from lsq to ode call; set up problem to call ode solver
global FIG_NUM;

kvar = kvar ./ p.scale;       %unscale parameters
p.kif = kvar(1);              %/s, forward CO2 hydration rate constant catalyzed by CA inside the coral tissue
p.kir = p.kif * (10^-p.pHe)./p.K1;              %uncatalyzed rate of HCO3 dehydration

CO2tot = [1 1 1] * p.cinit(1:3,1);           %total CO2ext conc at start of phase
Cinit(1:3,1) = p.cinit;
%calculate HCO3 concentrations assuming chemical and isotopic (apply binomial) equilibrium among species (this should be approximately correct in coral runs, validate) ;
Btot = p.DIC - CO2tot;                                        %HCO3- + CO32- concentration, call HCO3- and treat the same, making corrections to hyd/dehyd rates to account for it.
Cinit(4,1) = Btot * (1 - p.Taub)^3;                            %B61, row 4
Cinit(5,1) = Btot * (3 * (p.Taub) * (1-p.Taub)^2);            %B63, row 5
Cinit(6,1) = Btot * (3 * (p.Taub^2) * (1 - p.Taub));          %B65, row 6
Cinit(7,1) = Btot * (p.Taub^3);                                %B67, row 7

%approximate initial CO2 and HCO3 in the zoox as the same as external
Cinit(8:10,1)  = Cinit(1:3,1);
Cinit(11:14,1) = Cinit(4:7,1);

%use ode solver to deterimine time course of isotope exchange
options = odeset('RelTol', 1E-4, 'AbsTol', 1E-9);         %set ode options to obtain smooth (non-oscillating or diverging) solution
[t_ode, Cspecies] = ode15s(@coral_iCA_deriv, time, Cinit, options, p);
Cspecies = Cspecies';                  %ode gives one row per time point, transpose to return to normal structure: each row is the time series of a single Cspecies

  
if p.return == 0
CO2out = Cspecies(1:3,:);

elseif p.return == 1     
  CO2out = Cspecies;

  %write predicted data out to file, this is the fitted value
  outfile = 'isox_coral_iCA_fit.out';
  fout_id = fopen(outfile,'w');
  fprintf(fout_id,'time\t C45e\t C47e\t C49e\t B62e\t B64e\t B66e\t B68e\n');
  reps = size(Cspecies,2);
  for i = 1:reps
    fprintf(fout_id,'%6.2f\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\n',time(i),Cspecies(1:7,i));
  end
  fclose(fout_id);      
end

return
