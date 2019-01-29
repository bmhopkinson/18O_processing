function C13O2 = calibrate_C13O2(C13O2, p)
%calibrate and background subtract C13O2 data
%calibration assumes total CO2 is constant through the assay 

CO2tot = p.DIC ./ (1 + (p.K1./p.h) + ((p.K1 * p.K2)./(p.h^2)));      %calculate total 13CO2 from DIC and pH, see Zeebe and Wolf-Gladrow CO2 in SW Chp1 p4
fprintf(1,'CO2_total: %6.3E\n',CO2tot);

%BACKGROUND SUBTRACTION 
CO2_back = sum(C13O2(:,p.CYback_b:p.CYback_e),2);
CO2_back = CO2_back./(p.CYback_e - p.CYback_b + 1);   %calculate and subtract backgrounds for CO2_lab signals
C13O2 = C13O2 - (CO2_back * ones(1,p.cycles));

%CALIBRATION: use relative values for the dark portion 
sum_CO2 = [1; 1; 1] * sum(C13O2,1);		%add up all 13CO2 currents at each time point
C13O2 = CO2tot .* C13O2./sum_CO2;		%calibrated values




return