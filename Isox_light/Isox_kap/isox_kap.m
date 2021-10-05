function isox_kap(sample)
%processes isotope exchange data during photosynthesis
global FIG_NUM
FIG_NUM = 0;

file_param = strcat(sample,'.par');
file_data  = strcat(sample,'.cal');
file_ciup  = strcat(sample,'_Ciup.dat');

par  = load_params(file_param);
Ciup = load_Ciup(file_ciup);

%read in calibrated MS data
fdat = fopen(file_data,'r');
line = fgetl(fdat);                         %skip header line of file

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
    Ar(1,i) = A(6);                         %Argon, mass 40
    O2(1,i) = A(7);                         %Oxygen, mass 32    
    i = i+1;
end;
fclose(fdat);

tl = time(1,par.CYlight_b+1:par.CYlight_e-1);                     %trimming by 1 is necessary to match Ciup time
CO2_labl = CO2_lab(:,par.CYlight_b+1:par.CYlight_e-1);

%set up parameters for curve fitting
par.ce_init = CO2_labl(:,1);         %initial external CO2 values
par.O2_init = O2(par.CYlight_b+1);
%par.O2_init = 230E-9;
par.tl = tl;                %time used in interpolation of photo and Bup
par.Photo = Ciup.Photo;
par.Bup = Ciup.Bup;
par.flag = 0;           %return only 13CO2 data for fitting.


%call lsqcurvefit 
kcf_fit = kap_fit(tl,CO2_labl,par);
fprintf(1,'kcf_fit: %6.3E sd: %6.3E\n',kcf_fit(1,1),kcf_fit(1,2));
%kcf_fit  = 1.*par.kif;


%obtain full Cspecies and O2 time series based on fitted kcf
par.flag = 1;               %return full Cspecis and O2;
Cfit = kap_lsq(kcf_fit(1,1),tl,par);
t_ode = Cfit(1,:);

subplot(1,2,1)
plot(tl, CO2_labl,'o',t_ode, Cfit(2,:),'b',t_ode,Cfit(3,:),'g',t_ode,Cfit(4,:),'r'),title('extracelullar CO2');
subplot(1,2,2)
plot(t_ode,Cfit(9,:),'b',t_ode,Cfit(10,:),'g',t_ode,Cfit(11,:),'r'), title('intracellular CO2');

%figure(2)
%plot(tl,O2(par.CYlight_b+1:par.CYlight_e-1),t_ode, Cfit(16,:),'ob');

outfile ='isox_kap_fit.out';
fout = fopen(outfile, 'w');
fprintf(fout,'time\t C45_obs\t C47_obs\t C49_obs\t C45_fit\t C47_fit\t C49_fit\n');

for i= 1:length(tl)
    fprintf(fout,'%6.4f\t %6.4E\t %6.4E\t %6.4E\t %6.4E\t %6.4E\t %6.4E\n',tl(i),CO2_labl(:,i),Cfit(2:4,i));
end
fclose(fout);


return;
