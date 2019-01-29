function ydot = darkderiv(time, y, par)
%computes derivatives of C and H species for cells in the dark phase of O isotope exchange assay

par.kcr = par.kcf .* (10^-par.pHe)./par.K1;        %catalyzed rate of reverse rxn. in the dark assume pH ext = pH internal, so ratio of rates constant.

%correct hyd/dehyd rates to account for fraction of "h" pools that is actually HCO3-
par.kcr = par.kcr .* par.bicfrac;
par.kur = par.kur .* par.bicfrac;

%break out CO2 and HCO3 species
ce = y(1:3,1);                  %external CO2 species
he = y(4:7,1);                  %external HCO3 species
ci  = y(8:10,1);                 %internal CO2 species
hi  = y(11:14,1);                %internal HCO3 species

%define operator matrices         
H = zeros(3,4); 
H(1,1) = 1; H(1,2) = (1/3); H(2,2) = (2/3); H(2,3) = (2/3) ; H(3,3) = (1/3); H(3,4) = 1;                 %effect of H species on C species;
G =  eye(4,3);                          %effect of C species on H

%compute derivatives in mol/cm3/s
dce = -par.kuf .* ce + par.kur * H * he + (par.fc*par.cells)  .* (ci - ce);             %loss by hydration, gain by HCO3 dehydration, exchange w/ intracellular.
dhe =  par.kuf * G * ce - par.kur .* he + (par.fb*par.cells)  .* (hi - he);             %gain by hydration of CO2, loss by dehydration, exchange w/ intracellular. 
dci = -par.kcf .* ci + par.kcr * H * hi + (par.fc*par.cells/par.cyvol) .* (ce - ci);
dhi =  par.kcf * G * ci - par.kcr .* hi + (par.fb*par.cells/par.cyvol) .* (he - hi);

ydot(1:3,1)   = dce(1:3,1);
ydot(4:7,1)   = dhe(1:4,1);
ydot(8:10,1)  = dci(1:3,1);
ydot(11:14,1) = dhi(1:4,1);

return
