function Ciup = Ci_uptake(sample)
%calculates Ci uptake and photosynthesis as a function of time

file_param = strcat(sample,'.par');
file_data = strcat(sample,'.cal');

[par,init] = load_cipars(file_param);

fdat = fopen(file_data,'r');
line = fgetl(fdat);                         %skip header line of file

%initialize arrays
time    = zeros(1,par.cycles);
CO2_44  = zeros(1,par.cycles);
CO2_lab = zeros(3,par.cycles);
Ar      = zeros(1,par.cycles);
O2      = zeros(1,par.cycles);

%read in MS data from the rest of the file    
i = 1;
while ~feof(fdat)
    line = fgetl(fdat);
    A = sscanf(line, '%f %f %f %f %f %f %f');
    time(1,i) = A(1);                       %time in seconds
    CO2_44(1,i) = A(2);                     %mass 44 C12O16O16
    CO2_lab(1,i) = A(3);                    %mass 45 C13O16O16
    CO2_lab(2,i) = A(4);                    %mass 47 C13O18O16
    CO2_lab(3,i) = A(5);                    %mass 49 C13O18O18
    Ar(1,i) = A(6);                         % Argon, mass 40
    O2(1,i) = A(7);                         % Oxygen, mass 32    
    i = i + 1;
end;
fclose(fdat);

%filter O2 data during photosynthesis to smooth it

tp = time(par.CYlight_b:par.CYlight_e);
O2p = O2(par.CYlight_b:par.CYlight_e);
O2p= smooth(O2p,7);
O2p = O2p';             %smooth returns column vectors, change back to row
photo = ((O2p(3:end)-(O2p(1:end-2)))./(tp(3)-tp(1)))*par.evol;         %rate of photosynthesis in mol/s by finite differnce method

%truncate data to match size of photo
tp = tp(2:end-1);
O2p = O2p(2:end-1);
photos = photoFit(tp,photo);            %fit photosynthetic data with two lines representing a ramp up phase and steady-state phase.

figure(1)
subplot(3,2,1)
plot(time(par.CYlight_b:par.CYlight_e),O2(par.CYlight_b:par.CYlight_e),'og',tp,O2p), title('O2');
subplot(3,2,2)
plot(tp,photo,'or',tp,photos,'g'), title('photosynthesis')
photo = photos;

%determine net CO2 and HCO3 uptake
  %first track DIC and "HCO3-"
  DIC = ones(1,length(tp)).*par.DIC.*1E-9;            %units are mol/cm3
  delO2 = O2p - O2p(1);   %O2 produced during photosynthesis
  DIC = DIC - delO2;      %13C-DIC is lost to net carbon fixation 
  hp = 10^-par.pHe;       % H+ (mol/L; to match equilibrium constants)
  co2eq  = DIC ./ (1 + (par.K1./hp) + ((par.K1 * par.K2)./((hp)^2)));
  bicarb = DIC .* (par.K1 * hp)/(hp^2+par.K1*hp+par.K1*par.K2);   %bicarbonate
  carb = DIC .* (par.K1*par.K2)/(hp^2+par.K1*hp+par.K1*par.K2);   %carbonate
  be = bicarb + carb;             %term tracked in model is HCO3- + CO32- but called "HCO3-"

  %determine total 13CO2 and its rate of change
  tot13CO2 = sum(CO2_lab,1);              %total 13CO2 species
  tot13CO2p = smooth(tot13CO2(1,par.CYlight_b:par.CYlight_e),7);
  tot13CO2p = tot13CO2p';             %smooth returns column vectors, change back to row
  d13CO2p = ((tot13CO2p(3:end)-tot13CO2p(1:end-2))./(tp(3)-tp(1)));            %estimate derivative of tot13CO2

  %truncate tot13CO2p to match rate vectors
  tot13CO2p = tot13CO2p(2:end-1);

  %plot data
  subplot(3,2,3)
  plot(time(par.CYlight_b:par.CYlight_e),tot13CO2(par.CYlight_b:par.CYlight_e),'ob',tp,tot13CO2p, tp,co2eq,'m'),title('tot13CO2');
  subplot(3,2,4)
  plot(tp,d13CO2p,'og'),title('d13CO2/dt');

  %compute net CO2 uptake and net HCO3 uptake
  CO2up = par.evol.*(par.kur .* be - par.kuf .* tot13CO2p - d13CO2p);      %rate of net CO2 influx in mol/s
  subplot(3,2,5)
  plot(tp,CO2up,'or'), title('net CO2 uptake');

  Be_up = photo-CO2up;    %net HCO3- uptake in mol/s
%Be_up = smooth(Be_up);
%Be_up = Be_up';
subplot(3,2,6)
plot(tp,Be_up,'og'), title('net HCO3 uptake');


%write out data to file
outfile = strcat(sample,'_Ciup.dat');
fout = fopen(outfile,'w');

fprintf(fout,'time\t O2\t tot13CO2\t Photo\t CO2up\t HCO3up\n');
for i = 1:length(tp)
    fprintf(fout,'%6.4f\t %6.3E\t %6.3E\t %6.3E\t %6.3E\t %6.3E\n',tp(i),O2p(i),tot13CO2p(i),photo(i),CO2up(i),Be_up(i));
end
fclose(fout);

return;

    






