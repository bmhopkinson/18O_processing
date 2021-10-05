function ydot = coral_eCA_deriv(time,y,par)
%computes derivates of C and B species in 18O isotope exchnange due to
%coral

ce = y(1:3,1);          %extracellular CO2;
be = y(4:7,1);          %extracellular HCO3;
cs = y(8:10,1);         %surface CO2;
bs = y(11:14,1);        %surface HCO3-;
ci = y(15:17,1);        %intracellular CO2;
bi = y(18:21,1);        %intracellular HCO3-;

%define operator matrices         
H = zeros(3,4); 
H(1,1) = 1; H(1,2) = (1/3); H(2,2) = (2/3); H(2,3) = (2/3) ; H(3,3) = (1/3); H(3,4) = 1;   %effect of H species on C species;
G =  eye(4,3);                          %effect of C species on H

%correct dehydration rates for fraction of B that is CO32- rather than
%HCO3-
kur = par.kur .* par.bfrac;
ksr = par.ksr .* par.bfrac;
kcr = par.kcr .* par.bfrac;

if (par.ksfON == 0)
    ksf = 0;
    ksr = 0;
else
    ksf = par.ksf;
end

%compute derivatives in mol/cm3/s
dce = -par.kuf .* ce + kur * H * be + (par.Dc./par.DBL).*(par.Ac./par.Vl).*(cs - ce);
dbe =  par.kuf * G * ce - kur .* be + (par.Db./par.DBL).*(par.Ac./par.Vl).*(bs - be);
dcs = (-ksf./par.DBL) .* cs    + (ksr./par.DBL) * H * bs + (par.Dc./(par.DBL.^2)) .* (ce - cs) + (par.Pmc./par.DBL) .* (ci - cs);
dbs = ( ksf./par.DBL) * G * cs - (ksr./par.DBL) .* bs    + (par.Db./(par.DBL.^2)) .* (be - bs) + (par.Pmb./par.DBL) .* (bi - bs);
dci = -par.kcf .* ci + kcr * H * bi + (par.Pmc .* par.Ac ./par.Vc) .* (cs - ci);
dbi =  par.kcf * G * ci - kcr .* bi + (par.Pmb .* par.Ac ./par.Vc) .* (bs - bi);

ydot(1:3,1)   = dce(1:3,1);
ydot(4:7,1)   = dbe(1:4,1);
ydot(8:10,1)  = dcs(1:3,1);
ydot(11:14,1) = dbs(1:4,1);
ydot(15:17,1) = dci(1:3,1);
ydot(18:21,1) = dbi(1:4,1);

return

