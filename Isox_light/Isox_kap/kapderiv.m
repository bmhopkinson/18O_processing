function ydot = kapderiv(time, y, par)
%computes derivatives of C and H species for one compartment cells in the light phase of O isotope exchange assay

ce = y(1:3);
be = y(4:7);
ci = y(8:10);
bi = y(11:14);


%define operator matrices         
H = zeros(3,4); 
H(1,1) = 1; H(1,2) = (1/3); H(2,2) = (2/3); H(2,3) = (2/3) ; H(3,3) = (1/3); H(3,4) = 1;                 %effect of H species on C species;
G =  eye(4,3);                          %effect of C species on H

%transport should be proportional to relative abundance of CO2 masses;
ce_frac = ce ./ ([1 1 1] * ce);
be_frac = be ./ ([1 1 1 1] * be);
ci_frac  = ci  ./ ([1 1 1] * ci);

%rates of photosynthesis and HCO3 uptake at current time
%size(par.tl)
%size(par.Photo)
photo = interp1(par.tl, par.Photo, time, 'linear');
Bup   = interp1(par.tl, par.Bup, time, 'linear');

%photo = 0; %2.8E-17 .* par.cells;
%Bup = 0;%0.9E-17 .* par.cells;

ydot(1:3,1)   = -(par.kuf .* ce) + (par.kur .* H * be) + (par.fc.*par.cells./par.evol).*(ci - ce);
ydot(4:7,1)   =  (par.kuf .* G * ce) - (par.kur .* be) + (par.fb.*par.cells./par.evol).*(bi - be) - (Bup .* be_frac ./par.evol);
ydot(8:10,1)  = -(par.kcf_fit .* ci) + (par.kcr_fit .* H * bi) + (par.fc.*par.cells ./par.cyvol).*(ce - ci) - (photo./par.cyvol).*ci_frac;
ydot(11:14,1) =  (par.kcf_fit .* G *ci) - (par.kcr_fit .* bi) + (par.fb.*par.cells ./par.cyvol).* (be - bi) + (Bup./par.cyvol).* be_frac;
ydot(15)      =   photo;                  %O2 change


return
