function isox_eCA(infile)
%processes MIMS 18O isotope exchange data is the abscence of eCA inhibitor to provide eCA activity, requires paired experiment with eCA inhibitor
% to determine membrane permeabiliities and iCA activity
%v1.0

global FIG_NUM;
FIG_NUM = 0;

%READ IN PARAMETER AND DATA FILES
file_param = strcat(infile,'.par');                 %parameter file
data = strcat(infile,'.txt');                       %data file
par = load_params(file_param);                      %load parameter values 
[time, O2, Ar, CO2, C13O2] = load_data(data);     %read in MS dat
par.cycles = size(time,2);                          %number of data points

%calibrate 13CO3 data
C13O2 = calibrate_C13O2(C13O2,par);
tdat = time(1,par.CYdark_b:par.CYdark_e);       %time in which cells were present in the dark
CO2dat = C13O2(:,par.CYdark_b:par.CYdark_e);          %corresponding 13CO2 data   

%additional parameters for fitting, etc
par.infile = infile;
par.cinit = mean(C13O2(:,par.CYdark_b-3:par.CYdark_b),2);    %CO2 data at start of dark phase, average 3 points prior to cell addition
par.Taub = [0 0.5 1] * C13O2(:,par.CYuncat_e)./([1 1 1] * C13O2(:,par.CYuncat_e));      %18O atom fraction in HCO3.
par.plot = 0;                                      %plot the data? yes =1;
par.return = 0;                                    %return CO2 data only = 0; full Ci species =1;

%fit the data
kfit = fit_eCA(tdat, CO2dat,par);            %call routine that fits data collected w/ cells in the dark
se = kfit(:,1)-kfit(:,2);                          %mean - min;
fprintf(1,'ksf: %6.3E\t %6.3E\n', kfit(1,1), se(1));

%use fitted parameters to get a final predicted timecourse
par.plot = 1;                                      %plot data this time
par.return = 1;                                    %get full C species prediction
CO2pred = eCA_lsq(kfit, tdat, par);    %get predicted time course from fitted permeabilities and kcat

%plot results
FIG_NUM = FIG_NUM + 1;
figure(FIG_NUM)
subplot(2,2,1)
plot(tdat, CO2dat(1,:),'bo',tdat, CO2dat(2,:),'go', tdat, CO2dat(3,:),'ro', tdat, CO2pred(1,:),'b',  tdat, CO2pred(2,:),'g', tdat, CO2pred(3,:),'r'),title('CO2 bulk');
subplot(2,2,2)
plot(tdat, CO2pred(8,:),'b',tdat, CO2pred(9,:), 'g', tdat, CO2pred(10,:), 'r'),title('CO2 surface');
subplot(2,2,3)
plot(tdat, CO2pred(15,:),'b', tdat, CO2pred(16,:),'g',tdat, CO2pred(17,:),'r'),title('CO2 intracellular');
subplot(2,2,4)
plot(tdat, CO2pred(4,:),'b',tdat, CO2pred(5,:), 'g', tdat, CO2pred(6,:),'r', tdat, CO2pred(7,:), 'm', tdat, CO2pred(11,:),'b--',tdat, CO2pred(12,:),'g--',tdat, CO2pred(13,:),'r--',tdat,CO2pred(14,:),'m--'),title('HCO3- bulk and surface(dashed)');

FIG_NUM = FIG_NUM+1;
figure(FIG_NUM);
plot(tdat, CO2dat(1,:),'bo',tdat, CO2dat(2,:),'go', tdat, CO2dat(3,:),'ro', tdat, CO2pred(1,:),'b',  tdat, CO2pred(2,:),'g', tdat, CO2pred(3,:),'r'),title('CO2 bulk');
%append data needed for further analysis
fid = fopen(file_param,'a');
fprintf(fid,'ksf\t%e\n',kfit(1));
fclose(fid);

%write out "calibrated" CO2 data
CO2file = strcat(infile,'_eCA_CO2.out');
fCO2out = fopen(CO2file,'w');
fprintf(fCO2out,'time\t C45obs\t C47obs\t C49obs\n');
for i = 1:length(C13O2)
    fprintf(fCO2out,'%7.2f\t %6.3E\t %6.3E\t %6.3E\n',time(i),C13O2(:,i));
end
fclose(fCO2out);
    
return
