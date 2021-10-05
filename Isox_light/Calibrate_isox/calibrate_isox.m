function data = calibrate_isox(sample, FILETYPE)
%calibrate 13CO2 and O2 signals for analysis of light phase data
global FIG_NUM;
FIG_NUM = 0;
USE_AR = 0;

file_param = strcat(sample,'.par');
file_data = strcat(sample,'.txt');

[par,init] = load_params(file_param);

fdat = fopen(file_data,'r');

%read in MS data from the rest of the file    

if FILETYPE == 1
i = 1;
 while ~feof(fdat)
    line = fgetl(fdat);
    A = sscanf(line, '%i %f %f %f %f %f %f %f %f %f');
    time(1,i) = A(2);                       %time in seconds
    CO2_44(1,i) = A(4);                     %mass 44 C12O16O16
    CO2_lab(1,i) = A(5);                    %mass 45 C13O16O16
    CO2_lab(2,i) = A(6);                    %mass 47 C13O18O16
    CO2_lab(3,i) = A(7);                    %mass 49 C13O18O18
    Ar(1,i) = A(8);                         % Argon, mass 40
    O2(1,i) = A(9);                         % Oxygen, mass 32    
    i=i+1;
 end;
elseif FILETYPE == 2
  C = textscan(fdat, '%u %f %f %f %f %f %f %f %f');
  time = C{2}';
  O2 = C{4}';
  Ar = C{5}';
  CO2_44 = C{6}';
  CO2_lab = [C{7}';C{8}';C{9}'];
end

fclose(fdat);
par.cycles = length(O2);

%process O2 data
tscale = log((298.15 - (par.temp - 273.15))./par.temp);
O2sat = exp(5.80871 + 3.20291.*tscale + 4.17887.*(tscale^2) + 5.10006.*(tscale^3)+ (-9.86643E-2).*(tscale^4) + 3.80369.*(tscale^5)...
    + par.sal .* (-7.01577E-3 + (-7.70028E-3).*tscale + (-1.13864E-2).*(tscale^2) + (-9.51519E-3).*(tscale^3)) + (-2.75915E-7).*(par.sal^2));       %O2 saturation(umol/kg) from Pilson Intro to Chemistry of the Sea, corrected version of Garcia and Gordon LO 1992

if USE_AR == 0      %don't normalize O2 signal to Ar
O2signal_sat = sum(O2(par.CYuncat_b:par.CYuncat_e),2)./length(O2(par.CYuncat_b:par.CYuncat_e));
cfO2 = O2sat./O2signal_sat;         %calibration factor for O2 current to O2 concentration
O2 = O2 .* cfO2 .* 1E-9;            %O2 in mol/cm3

elseif USE_AR == 1      %normalize O2 signal to Ar
O2Ar = O2./Ar;
O2Ar_sat = sum(O2Ar(par.CYuncat_b:par.CYuncat_e),2)./length(O2Ar(par.CYuncat_b:par.CYuncat_e));
cfO2 = O2sat./O2Ar_sat;                  %calibration factor for O2/Ar -> O2
O2 = (O2Ar.*cfO2).*1E-9;                 %O2 in mol/cm3

end

%plot O2 data
FIG_NUM = FIG_NUM +1;
figure(FIG_NUM);
plot(time,O2,'b'),title('O2 (mol/cm3)');


%process 13CO2 data
%calculate calibration factors for 13CO2 species assuming CO2 in uncatlyzed
%and dark phase in chemical equilibrium with HCO3

%equilbrium 13CO2 concentration
h = 10^-par.pHe;            %H+ conc
CO2eq = (par.DIC*1E-6) ./ (1 + (par.K1./h) + ((par.K1 * par.K2)./(h^2)));      %calculate total 13CO2 from DIC and pH, see Zeebe and Wolf-Gladrow CO2 in SW Chp1 p4
CO2eq = CO2eq./(1000);                           %convert to units of mol/cm3

CO2back = sum(CO2_lab(:,par.CYback_b:par.CYback_e),2)./(par.CYback_e - par.CYback_b + 1);   %calculate and subtract backgrounds for CO2_lab signals
CO2_lab = CO2_lab - (CO2back * ones(1,par.cycles));

CO2uncat = CO2_lab(:,par.CYuncat_b:par.CYuncat_e);
us = size(CO2uncat,2);
CO2dark = CO2_lab(:,par.CYdark_b+10:par.CYdark_e);  % CYDARK_b+10 to skip rapid phase
CO2 = cat(2,CO2uncat,CO2dark);
cfs = calibration(CO2, CO2eq, par.enrich, us);
cfCO2 = diag(cfs);          %matrix with calibration factors on the diagonal, right now all the CFs are the same, its the best option right now.
CO2calib = cfCO2 * CO2_lab;

%if there is a light phase, rescale total 13CO2 so that it is in equilbrium at the start of the light phase
if ~(par.CYlight_b == -999)
tot13CO2 = sum(CO2calib,1);              %total 13CO2 species
rescale = CO2eq./mean(tot13CO2(par.CYlight_b-3:par.CYlight_b));
CO2calib = rescale .* CO2calib;
end

%plot calibrated 13CO2 data
FIG_NUM = FIG_NUM +1;
figure(FIG_NUM);
plot(time,CO2calib),title('13CO2');

%write out CO2 concentration data
calib_out = strcat(sample,'.cal');           %.cal extension for calibrated 13CO2 + additional signals output file.
fnt = fopen(calib_out,'w');

fprintf(fnt,'time\t C44\t C45\t C47\t C49\t Ar\t O2\n');
steps = length(time);
for i = 1:steps
    fprintf(fnt,'%6.2f\t %6.3E\t %6.3E\t %6.3E\t %6.3E\t %6.3E\t %6.3E\n',time(i),CO2_44(1,i),CO2calib(:,i),Ar(1,i),O2(1,i));
end;
fclose(fnt);
return


