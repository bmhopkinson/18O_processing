function isox_uncat(infile, FILETYPE)
%processes MIMS 18O isotope exchange data to find uncatalyzed CO2
%hydration/dehydration rates 
%v1.0

global FIG_NUM;
FIG_NUM = 0;
file_param = strcat(infile,'.par');
file_data  = strcat(infile,'.txt');

%READ IN DATA FILE: PARAMETERS AND MS DATA
par = load_params(file_param);
[time, O2, Ar, CO2, C13O2] = load_data(file_data,FILETYPE);     %read in MS data 

%compute several additional parameters 
par.cycles = size(time,2);
par.h = 10^-par.pHe;            %H+ conc
par.DIC = par.DIC./(1000);     %convert to units of mol/cm3

%calculate predicted kinetic constants, terminology of Zeebe and
%Wolf-Gladrow CO2 in seawater book
kp1 = exp(1246.98 - 6.19E4 ./ par.temp - 183 .* log(par.temp));         % Johnson 1982
kp4 = 4.70E7 .* exp(-23.2 ./ (8.314E-3 .* par.temp));                   % Johnson 1982

kuf_exp = kp1 + kp4 .* (par.Kw ./(10^-par.pHe));                        % effective kuf expected  = kp1 + kp4 * [OH-]
fprintf(1,'expected kuf: %6.4E\n',kuf_exp);

%calibrate 13CO3 data
C13O2 = calibrate_C13O2(C13O2,par);  

par.cinit = mean(C13O2(:,(par.CYuncat_b:par.CYuncat_b+3)),2);       %initial values for ode

%ANALYZE UNCATALZED PHASE
%calculate uncatalyzed interconversion rate 
kuf = fit_uncat(time(1,par.CYuncat_b:par.CYuncat_e), C13O2(:,par.CYuncat_b:par.CYuncat_e), par);   
kur = kuf * (10^-par.pHe)./par.K1;
fprintf(1,'kuf_fit: %6.4E\n',kuf);      %write out result to the screen

%append kuf value to parameter file
fpar = fopen(file_param,'a');            %open file
fprintf(fpar,'kuf\t%6.4E\n',kuf);
fclose(fpar);

return
