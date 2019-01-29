function isox_dark(infile, FILETYPE)
%processes MIMS 18O isotope exchange data to provide cellular transfer coefficients and intracellular CA activity
%v1.0

global FIG_NUM;
FIG_NUM = 0;

%READ IN PARAMETER AND DATA FILES
file_param = strcat(infile,'.par');          %parameter file
data = strcat(infile,'.txt');            %data file
par = load_params(file_param);           %load parameter values (initial ci concentrations not used here).
[time, O2, Ar, CO2, C13O2] = load_data(data,FILETYPE);     %read in MS data 

%compute several additional parameters 
par.cycles = size(time,2);
par.h = 10^-par.pHe;            %H+ conc
par.DIC = par.DIC./(1000);                                 %convert from M to mol/cm3

%calibrate 13CO3 data
C13O2 = calibrate_C13O2(C13O2,par);

%ANALYZE DARK PHASE
time_cd = time(1,par.CYdark_b:par.CYdark_e);                %time in which cells were present in the dark
CO2_lab_cd = C13O2(:,par.CYdark_b:par.CYdark_e);          %corresponding 13CO2 data

%additional parameters for fitting
par.cinit = CO2_lab_cd(1:3,1);                  %CO2 data at start of dark phase
par.Taub = [0 0.5 1] * C13O2(:,par.CYuncat_e)./([1 1 1] * C13O2(:,par.CYuncat_e));      %18O atom fraction in HCO3.
par.plot = 0;               %plot the data? yes =1;
par.return = 0;               %return CO2 data only = 0; full Ci species =1;
par.scale = [0.001 1E8 1E10];     %scale parameters to be fit to help optimization algorithim

kfit = fit_dark(time_cd, CO2_lab_cd,par);       %call routine which fits data collected w/ cells in the dark

%print out fitted parameters and confidence intervals
fprintf(1,'kcf (/s): %e\t sd: %e\n',kfit(1,1),kfit(1,2));
fprintf(1,'fc (cm3/s): %e\t sd: %e\n',kfit(2,1),kfit(2,2));
fprintf(1,'fb (cm3/s): %e\t sd: %e\n',kfit(3,1),kfit(3,2));


%use fitted parameters to get a final predicted timecourse
par.plot = 1;                                    %plot data this time
par.return = 1;                                    %get full C species prediction

CO2pred_cd = dark_lsq(kfit(:,1)'.*par.scale, time_cd, par);      %get predicted time course from fitted permeabilities and kcat
FIG_NUM = FIG_NUM + 1;
figure(FIG_NUM)
plot(time_cd, CO2_lab_cd(1,:),'bo',time_cd, CO2_lab_cd(2,:),'go', time_cd, CO2_lab_cd(3,:),'ro', time_cd, CO2pred_cd(1,:),'b',  time_cd, CO2pred_cd(2,:),'g', time_cd, CO2pred_cd(3,:),'r'),title('Dark Phase Fit vs Data');


%append data needed for further analysis
fid = fopen(file_param,'a');
fprintf(fid,'kcf\t%e\n',kfit(1));
fprintf(fid,'fc\t%e\n',kfit(2));
fprintf(fid,'fb\t%e\n',kfit(3));
fprintf(fid,'HCO3_ext\t%e\t %e\t %e\t %e\n',CO2pred_cd(4:7,end));
fprintf(fid,'CO2_in\t %e\t %e\t %e\n',CO2pred_cd(8:10,end));
fprintf(fid,'HCO3_in\t %e\t %e\t %e\t %e\n',CO2pred_cd(11:14,end));

%write fitted values to outfile
outfile = strcat(infile,'_passive_params.out');
fout = fopen(outfile, 'w');
fprintf(fout, 'parameter\tmean\tstd error\n');
fprintf(fout, 'kcf\t%e\t%e\n',kfit(1,1),kfit(1,2));
fprintf(fout, 'fc\t%e\t%e\n',kfit(2,1), kfit(2,2));
fprintf(fout, 'fb\t%e\t%e\n',kfit(3,1), kfit(3,2));
fclose(fout);


return
