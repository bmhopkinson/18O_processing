function ydot = coral_iCA_deriv(time,y,p)
%computes derivates of C and B species in 18O isotope exchnange due to
%coral

ce = y(1:3,1);          %extracellular CO2;
be = y(4:7,1);          %extracellular HCO3;
ci = y(8:10,1);         %intracellular CO2;
bi = y(11:14,1);        %intracellular HCO3-;

%define operator matrices         
H = zeros(3,4); 
H(1,1) = 1; H(1,2) = (1/3); H(2,2) = (2/3); H(2,3) = (2/3) ; H(3,3) = (1/3); H(3,4) = 1;                 %effect of H species on C species;
G =  eye(4,3);                          %effect of C species on H

kif_s = p.kif .* (p.Vc .* p.Vadd)./(p.Vhom .* p.Ve);
kir_s = p.bfrac .* kif_s .* p.h./p.K1; 

%compute derivatives in mol/cm3/s
dce = -(p.kuf + p.ksf_s + kif_s) .* ce + (p.kur + p.ksr_s + kir_s) * H * be + (p.fc.* p.cells).*(ci - ce);
dbe =  (p.kuf + p.ksf_s + kif_s) * G * ce - (p.kur + p.ksr_s + kir_s) .* be + (p.fb.* p.cells).*(bi - be);
dci = -p.kcf .* ci + p.kcr * H * bi + (p.fc ./p.Vcell) .* (ce - ci);
dbi =  p.kcf * G * ci - p.kcr .* bi + (p.fb ./p.Vcell) .* (be - bi);

ydot(1:3,1)   = dce(1:3,1);
ydot(4:7,1)   = dbe(1:4,1);
ydot(8:10,1)  = dci(1:3,1);
ydot(11:14,1) = dbi(1:4,1);

return

