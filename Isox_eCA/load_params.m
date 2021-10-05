function par = load_params(file)
%opens parameter file and loads in data 

fid = fopen(file,'r');
if fid==-1
   error(['File ' file ' not found or permission denied.']);
end

ct = 1;
while ~feof(fid)
    raw{ct}=fgetl(fid);
    ct=ct+1;
end
fclose(fid);

%define all required parameters to be read from file in a structure and set as empty. later will
%check if any are left empty 
par = struct('DIC',[],'evol',[],'pHe',[],'temp',[],'cells',[],'Reff',[],'cyvol',[],'cySA',[],'enrich',[],...
    'CYback_b',[],'CYback_e',[],'CYuncat_b',[],'CYuncat_e',[],'CYdark_b',[],'CYdark_e',[],'CYlight_b',[],...
    'CYlight_e',[],'kuf',[],'kcf',[],'fc',[],'fb',[]);

for i = 1:length(raw)
    [id] = sscanf(raw{i},'%s[\t]%*s');
    
    switch id
        
        case 'DIC'
            par.DIC = (sscanf(raw{i},'%*s %f'))*1E-6;       %DIC in M
        case 'evol'
            par.evol = sscanf(raw{i},'%*s %f');
        case 'pHe'
            par.pHe = sscanf(raw{i},'%*s %f');
        case 'temp'
            par.temp = 273.15 + sscanf(raw{i},'%*s %f');        %temp in K
        case 'cells/mL'
            par.cells = sscanf(raw{i},'%*s %e');
        case 'Rshape'
            par.Reff = sscanf(raw{i},'%*s %f');
            par.Reff = par.Reff * (100/1E6);        %shape radius converted from um to cm
        case 'cell_vol'
            par.cyvol = sscanf(raw{i},'%*s %e');
        case 'cell_SA'
            par.cySA = sscanf(raw{i},'%*s %e');
        case 'enrich_factor'
            par.enrich = sscanf(raw{i},'%*s %f');
        case 'background_begin'
            par.CYback_b = sscanf(raw{i},'%*s %d');
        case 'background_end'
            par.CYback_e = sscanf(raw{i},'%*s %d');
        case 'uncat_begin'
            par.CYuncat_b = sscanf(raw{i},'%*s %d');
        case 'uncat_end'
            par.CYuncat_e = sscanf(raw{i},'%*s %d');
        case 'dark_begin'
            par.CYdark_b = sscanf(raw{i},'%*s %d');
        case 'dark_end'
            par.CYdark_e = sscanf(raw{i},'%*s %d');
        case 'light_begin'
            par.CYlight_b = sscanf(raw{i},'%*s %d');
        case 'light_end'
            par.CYlight_e = sscanf(raw{i},'%*s %d');
        case 'kuf'
            par.kuf = sscanf(raw{i},'%*s %e');
        case 'kcf'
            par.kcf = sscanf(raw{i},'%*s %e');
        case 'fc'
            par.fc = sscanf(raw{i},'%*s %e');
        case 'fb'
            par.fb =  sscanf(raw{i},'%*s %e');
    end
end

%check to make sure all the fields were filled
fn = fieldnames(par);
for i = 1:length(fn)
    if isempty(par.(fn{i}))
        error('Error in load_params, missing parameter: %s\n',fn{i})
    end
end

%additional parameters
par.TC = par.temp - 273.15;
par.h = 10^-par.pHe;            %H+ conc
par.DIC = par.DIC./(1000);                                 %convert from M to mol/cm3
par.sal = 35;                   %salinity assumed to be 35, doesn't have large effect on K1.                   
par.K1 = 10^-((3633.86./par.temp) - 61.2172 + 9.6777*log(par.temp) - 0.011555 .* par.sal + 0.0001152 .* (par.sal^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
par.K2 = 10^((-471.78/par.temp) - 25.929 + 3.16967 * log(par.temp) + 0.01781 * par.sal - 0.0001122*par.sal^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
par.Kw = exp(148.96502 - (13847.26 ./ par.temp) - 23.6521 .* log(par.temp) + (par.sal^.5).*((118.67 ./ par.temp) - 5.977 + 1.0495 .* log(par.temp)) - 0.01615 .* par.sal); %ion product of water, CO2 methods DOE 1994

%modify parameters for use in model
par.cyvol = par.cells * par.cyvol;                %one-box model for dark phase, use total cell volume/mL;
par.cySA = par.cells * par.cySA;
par.kur = par.kuf .* par.h./par.K1;              %uncatalyzed rate of HCO3 dehydrationpar.
par.kcr = par.kcf .* par.h./par.K1;     %catalyzed rate of reverse rxn. in the dark assume pH ext = pH internal, so ratio of rates constant.

%calculate mass transfer coefficients and surface volume
par.ls   = 1E-5; %radial length of surface layer
DCF = diff_coef(par.TC,par.sal, 1);        %diffusion coefficients (in cm2/s) at current temp, salinity, final variable is pressure: assume 1 atm.
par.Dc  = DCF(3);     %diffusivity of CO2 in cm2/s
par.Db  = DCF(8);       %diffusivity of HCO3- in cm2/s
par.fcBS = 4*pi * par.Dc * (par.Reff + par.ls); %parameter for diffusive CO2 flux to cell surface layer
par.fbBS = 4*pi * par.Db * (par.Reff + par.ls); %parameter for diffusive HCO3- flux to cell surface layer
par.fcSM = 1/((1/par.fc)-(1/par.fcBS));      %mass transfer coeff for CO2 from boundary to surface layer
par.fbSM = 1/((1/par.fb)-(1/par.fbBS));      %HCO3- membrane permeability flux parameter (Pb*A)
par.svol = (4/3)*pi*((par.Reff+par.ls)^3 - par.Reff^3); % effective surface volume in cm3

par;

return