function isox_coral_eCA(infile)
%process MIMS 18O exchange data to provide eCA activity; also fits internal
%CA activity and membrane permeabilites from data obtained in the presence
%of DBAZ, but these parameters are not well constrained and should be used
%only to account for the effects of iCA activity earlier in the run (prior
%to addition of DBAZ). 
global FIG_NUM;
FIG_NUM = 0;


%read in parameter and data files
file_param = strcat(infile,'.par');             %parameter file
data = strcat(infile,'.txt');                   %data file
par = load_params(file_param);                  %load parameter values (initial ci concentrations not used here).
[time, O2, Ar, CO2, C13O2] = load_data(data);   %read in MS data 

%define several additional parameters 
par.infile = infile;
par.cycles = size(time,2);        
par.h = 10^-par.pHe;            %H+ conc
par.DIC = par.DIC./(1000);                                 %convert from M to mol/cm3

%calibrate 13CO3 data
C13O2 = calibrate_C13O2(C13O2,par);
dlmwrite('C13O2_calib.txt',[time' C13O2'],'\t');

%ANALYZE DARK PHASE
time_CA = time(1,par.CYeCA_b:par.CYiCA_e-3);          %eCA phase: coral in chamber, no inhibitors
CO2_CA = C13O2(:,par.CYeCA_b:par.CYiCA_e-3);          %corresponding 13CO2 data

%additional parameters for fitting
par.brk = par.CYeCA_e - par.CYeCA_b +1;                %"break" in data segments based on addition of DBAZ 
par.ci1 = sum(CO2_CA(1:3,1:3),2)./3;                    %initial CO2 data at start of dark phase
par.ci2 = sum(CO2_CA(1:3,par.brk+1:par.brk+3),2)./3;    %initial CO2 data after DBAZ addition
par.Taub = [0 0.5 1] * par.ci1./([1 1 1] * par.ci1);      %18O atom fraction in HCO3, need to validate that CO2 18O content is ~ HCO3- 18O content in coral runs.
par.teCAend = time(par.CYeCA_e);
par.plot = 0;               %plot the data? yes =1;
par.return = 0;               %return CO2 data only = 0; full Ci species =1;
par.scale = [1 0.1 1E2 1E2];     %scale parameters to be fit to help optimization algorithim

kfit = fit_coral_eCA(time_CA, CO2_CA,par);       %call routine which fits data collected w/ cells in the dark

%print out fitted parameters and confidence intervals
fprintf(1,'ksf (cm/s): %e\t sd: %e\n',kfit(1,1),kfit(1,2));
fprintf(1,'kcf (/s): %e\t sd: %e\n',kfit(2,1),kfit(2,2));
fprintf(1,'Pc (cm/s): %e\t sd: %e\n',kfit(3,1),kfit(3,2));
fprintf(1,'Pb (cm/s): %e\t sd: %e\n',kfit(4,1),kfit(4,2));

%use fitted parameters to get a final predicted timecourse
par.plot = 1;                                    %plot data this time
par.return = 1;                                    %get full C species prediction

CO2pred = coral_eCA_lsq(kfit(:,1)'.*par.scale, time_CA, par);      %get predicted time course from fitted permeabilities and kcat
FIG_NUM = FIG_NUM + 1;
figure(FIG_NUM)
plot(time_CA, CO2_CA(1,:),'bo',time_CA, CO2_CA(2,:),'go', time_CA, CO2_CA(3,:),'ro', time_CA, CO2pred(1,:),'b',  time_CA, CO2pred(2,:),'g', time_CA, CO2pred(3,:),'r'),title('eCA/iCA Fit vs Data');

return
