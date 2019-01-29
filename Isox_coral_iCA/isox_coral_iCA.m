function isox_coral_iCA(infile)
%process MIMS 18O exchange data of coral homogenate to determine internal
%CA activity of coral tissue (accounts for zoox contamination and eCA
%activity)
global FIG_NUM;
FIG_NUM = 0;

%read in parameter and data files
file_param = strcat(infile,'.par');             %parameter file
data = strcat(infile,'.txt');                   %data file
[index,time, O2, Ar, CO2, C13O2] = load_data(data);   %read in MS data
p = load_params(file_param, index);                  %load parameter values (initial ci concentrations not used here).

%compute several additional parameters 
p.cycles = size(time,2);

%calibrate 13CO3 data
C13O2 = calibrate_C13O2(C13O2,p);
dlmwrite('C13O2_calib.txt',[time' C13O2'],'\t');  %output calibrated data

%ANALYZE DARK PHASE
time_hom = time(1,p.CYhom_b:p.CYhom_e);          %iCA phase: coral+eCA inhibitor
CO2_hom = C13O2(:,p.CYhom_b:p.CYhom_e);          %corresponding 13CO2 data

%additional parameters for fitting
p.cinit = mean(CO2_hom(1:3,1:2),2);                  %CO2 data at start of dark phase
p.Taub = [0 0.5 1] * p.cinit./([1 1 1] * p.cinit);      %18O atom fraction in HCO3, need to validate that CO2 18O content is ~ HCO3- 18O content in coral runs.
p.plot = 0;               %plot the data? yes =1;
p.return = 0;             %return CO2 data only = 0; full Ci species =1;
p.scale = 1;        %scale parameters to be fit to help optimization algorithim

kfit = fit_coral_iCA(time_hom, CO2_hom,p);       %call routine which fits data collected w/ cells in the dark

%print out fitted parameters and confidence intervals
fprintf(1,'kif (/s): %e\t sd: %e\n',kfit(1,1),kfit(1,2));

%use fitted parameters to get a final predicted timecourse
p.plot = 1;           %plot data this time
p.return = 1;         %get full C species prediction

CO2pred = coral_iCA_lsq(kfit(1,1).*p.scale, time_hom, p);      %get predicted time course from fitted permeabilities and kcat
FIG_NUM = FIG_NUM + 1;
figure(FIG_NUM)
plot(time_hom, CO2_hom(1,:),'bo',time_hom, CO2_hom(2,:),'go', time_hom, CO2_hom(3,:),'ro', time_hom, CO2pred(1,:),'b',  time_hom, CO2pred(2,:),'g', time_hom, CO2pred(3,:),'r'),title('Homogenate Phase Fit vs Data');

return
