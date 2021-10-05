function ydot = eCAderiv(time, y, par)
%computes derivatives of C and H species for cells in the dark phase of O isotope exchange assay

par.ksr = par.ksf .* (10^-par.pHe)./par.K1;     %catalyzed rate of reverse rxn. in the dark assume pH ext = pH internal, so ratio of rates constant.

%correct hyd/dehyd rates to account for fraction of "h" pool that is actually HCO3-
fracbic = (1/ (1+ par.K2./(10^-par.pHe)));           % fraction of "h" pool that is HCO3-
par.kur = par.kur .* fracbic;
par.ksr = par.ksr .* fracbic;
par.kcr = par.kcr .* fracbic;


%break out CO2 and HCO3 species
ce = y(1:3,1);                   %external CO2 species
be = y(4:7,1);                   %external HCO3 species
cs = y(8:10,1);                  %surface CO2 species
bs = y(11:14,1);                 %surface HCO3- species
ci = y(15:17,1);                 %intracellular CO2 species
bi = y(18:21,1);                %intracellular HCO3 species

%define operator matrices         
H = zeros(3,4); 
H(1,1) = 1; H(1,2) = (1/3); H(2,2) = (2/3); H(2,3) = (2/3) ; H(3,3) = (1/3); H(3,4) = 1;                 %effect of H species on C species;
G =  eye(4,3);                          %effect of C species on H

%compute derivatives in mol/cm3/s
dce = -par.kuf .* ce + par.kur * H * be + (par.fcBS*par.cells) .* (cs - ce);             %loss by hydration, gain by HCO3 dehydration, exchange w/ intracellular.
dbe =  par.kuf * G * ce - par.kur .* be + (par.fbBS*par.cells) .* (bs - be);             %gain by hydration of CO2, loss by dehydration, exchange w/ intracellular. 
dcs = (par.cells./par.svol)*(-par.ksf .* cs + par.ksr * H * bs  + (par.fcBS) .* (ce - cs) + (par.fcSM) .* (ci - cs));
dbs = (par.cells./par.svol)*( par.ksf .* G * cs - par.ksr .* bs + (par.fbBS) .* (be - bs) + (par.fbSM) .* (bi - bs));
dci = -par.kcf .* ci + par.kcr * H * bi + (par.fcSM*par.cells/par.cyvol) .* (cs - ci);
dbi =  par.kcf * G * ci - par.kcr .* bi + (par.fbSM*par.cells/par.cyvol) .* (bs - bi);

ydot(1:3,1)   = dce(1:3,1);
ydot(4:7,1)   = dbe(1:4,1);
ydot(8:10,1)  = dcs(1:3,1);
ydot(11:14,1) = dbs(1:4,1);
ydot(15:17,1) = dci(1:3,1);
ydot(18:21,1) = dbi(1:4,1);

return
