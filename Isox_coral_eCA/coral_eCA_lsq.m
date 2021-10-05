function CO2out = coral_eCA_lsq(kvar, time, par)
%interface from lsq to ode call; set up problem to call ode solver
global FIG_NUM;

kvar = kvar ./ par.scale;       %unscale parameters
par.ksf = kvar(1);              %/s, forward CO2 hydration rate constant catalyzed by CA inside the cell, initial guess
par.kcf = kvar(2);              %/s, rate constant for diffusive CO2 flux
par.Pmc = kvar(3);              %CO2 permeability 
par.Pmb = kvar(4);              %HCO3- permeability assumed to be zero

par.kcr = par.kcf * (10^-par.pHe)./par.K1;              %intracellular rate of HCO3 dehydration
par.ksr = par.ksf * (10^-par.pHe)./par.K1;              %surface eCA catalyzed rate of HCO3- dehydration

CO2tot = [1 1 1] * par.ci1(1:3,1);           %total CO2ext conc at start of phase
Cinit(1:3,1) = par.ci1;
%calculate HCO3 concentrations assuming chemical and isotopic (apply binomial) equilibrium among species (this should be approximately correct in coral runs, validate) ;
Btot = par.DIC - CO2tot;                                        %HCO3- + CO32- concentration, call HCO3- and treat the same, making corrections to hyd/dehyd rates to account for it.
Cinit(4,1) = Btot * (1 - par.Taub)^3;                           %B61, row 4
Cinit(5,1) = Btot * (3 * (par.Taub) * (1-par.Taub)^2);          %B63, row 5
Cinit(6,1) = Btot * (3 * (par.Taub^2) * (1 - par.Taub));        %B65, row 6
Cinit(7,1) = Btot * (par.Taub^3);                               %B67, row 7

%approximate initial CO2 and HCO3 in the boundary layer and inside the coral as equivalant to external
Cinit(8:10,1)  = Cinit(1:3,1);
Cinit(11:14,1) = Cinit(4:7,1);
Cinit(15:17,1) = Cinit(1:3,1);
Cinit(18:21,1) = Cinit(4:7,1);

%eCA portion use ode solver to deterimine time course of isotope exchange
par.ksfON = 1;
t1 = time(1,1:par.brk);
options = odeset('RelTol', 1E-4, 'AbsTol', 1E-9);         %set ode options to obtain smooth (non-oscillating or diverging) solution
[t_ode, Cspecies] = ode15s(@coral_eCA_deriv, t1, Cinit, options, par);
Cspecies = Cspecies';                  %ode gives one row per time point, transpose to return to normal structure: each row is the time series of a single Cspecies
%iCA portion
par.ksfON = 0;
t2 = time(1,par.brk+1:end);
Cinit = [par.ci2; Cspecies(4:end,end)];
options = odeset('RelTol', 1E-4, 'AbsTol', 1E-9);         %set ode options to obtain smooth (non-oscillating or diverging) solution
[t_ode, Cspecies2] = ode15s(@coral_eCA_deriv, t2, Cinit, options, par);
Cspecies = [Cspecies Cspecies2'];

if par.return == 0
CO2out = Cspecies(1:3,:);

elseif par.return == 1     
CO2out = Cspecies;
  
  %write predicted data out to file, this is the fitted value
  outfile = strcat(par.infile,'_eCAfit.txt');
  fout_id = fopen(outfile,'w');
  fprintf(fout_id,'time\t C45e\t C47e\t C49e\t B62e\t B64e\t B66e\t B68e\t C45s\t C47s\t C49s\t B62s\t B64s\t B66s\t B68s\t C45i\t C47i\t C49i\t B62i\t B64i\t B66i\t B68i\n');
  reps = size(Cspecies,2);
  for i = 1:reps
    fprintf(fout_id,'%6.2f\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\t %5.3E\n',time(i),Cspecies(:,i));
  end
  fclose(fout_id);     
  
  if par.plot == 1
      FIG_NUM = FIG_NUM + 1;
      figure(FIG_NUM)
      subplot(2,3,1)
      plot(time, Cspecies(1,:),'b',time, Cspecies(2,:), 'g',time, Cspecies(3,:), 'r'),title('bulk CO2');
      subplot(2,3,2)
      plot(time, Cspecies(8,:),'b',time, Cspecies(9,:), 'g',time, Cspecies(10,:), 'r'),title('surface CO2');
      subplot(2,3,3)
      plot(time, Cspecies(15,:),'b',time, Cspecies(16,:), 'g',time, Cspecies(17,:), 'r'),title('intracellular CO2');
      subplot(2,3,4)
      plot(time, Cspecies(4,:),'b',time, Cspecies(5,:), 'g',time, Cspecies(6,:), 'r', time, Cspecies(7,:), 'm'),title('bulk HCO3-');
      subplot(2,3,5)
      plot(time, Cspecies(11,:),'b',time, Cspecies(12,:), 'g',time, Cspecies(13,:), 'r', time, Cspecies(14,:), 'm'),title('surface HCO3-');
      subplot(2,3,6)
      plot(time, Cspecies(18,:),'b',time, Cspecies(19,:), 'g',time, Cspecies(20,:), 'r', time, Cspecies(21,:), 'm'),title('intracellular HCO3-');
  end
  
      
  
end

return
