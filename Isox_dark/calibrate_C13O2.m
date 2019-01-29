function C13O2 = calibrate_C13O2(C13O2, par)
%calibrate and background subtract C13O2 data

CO2tot = par.DIC ./ (1 + (par.K1./par.h) + ((par.K1 * par.K2)./(par.h^2)));      %calculate total 13CO2 from DIC and pH, see Zeebe and Wolf-Gladrow CO2 in SW Chp1 p4
fprintf(1,'CO2_total: %6.3E\n',CO2tot);
fprintf(1,'cells %4E, cell_area %6.3E, cell_vol %6.3E\n',par.cells, par.cySA ,par.cyvol);

%BACKGROUND SUBTRACTION 
CO2_back = sum(C13O2(:,par.CYback_b:par.CYback_e),2)./(par.CYback_e - par.CYback_b + 1);   %calculate and subtract backgrounds for CO2_lab signals
C13O2 = C13O2 - (CO2_back * ones(1,par.cycles));

%CALIBRATION: use relative values for the dark portion 
sum_CO2 = [1; 1; 1] * sum(C13O2(:,par.CYuncat_b:par.CYdark_e),1);		%add up all 13CO2 currents at each time point
C13O2(:,par.CYuncat_b:par.CYdark_e) = CO2tot .* (C13O2(:,par.CYuncat_b:par.CYdark_e)./sum_CO2);		%calibrated values

return