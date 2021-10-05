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
par = struct('DIC',[],'Vl',[],'pHe',[],'temp',[],'Ac',[],'tx',[],'DBL',[],...
    'CYback_b',[],'CYback_e',[],'CYuncat_b',[],'CYuncat_e',[],'CYeCA_b',[],'CYeCA_e',[],'CYiCA_b',[],...
    'CYiCA_e',[],'kuf',[]);

for i = 1:length(raw)
    [id] = sscanf(raw{i},'%s[\t]%*s');
    
    switch id
        
        case 'DIC'
            par.DIC = (sscanf(raw{i},'%*s %f'))*1E-6;       %DIC in M
        case 'evol'
            par.Vl = sscanf(raw{i},'%*s %f');
        case 'pHe'
            par.pHe = sscanf(raw{i},'%*s %f');
        case 'temp'
            par.temp = 273.15 + sscanf(raw{i},'%*s %f');        %temp in K
        case 'SA'
            par.Ac = sscanf(raw{i},'%*s %e');
        case 'tissue_thickness'
            par.tx = sscanf(raw{i},'%*s %e');
        case 'DBL'
            par.DBL = sscanf(raw{i},'%*s %e');
        case 'background_begin'
            par.CYback_b = sscanf(raw{i},'%*s %d');
        case 'background_end'
            par.CYback_e = sscanf(raw{i},'%*s %d');
        case 'uncat_begin'
            par.CYuncat_b = sscanf(raw{i},'%*s %d');
        case 'uncat_end'
            par.CYuncat_e = sscanf(raw{i},'%*s %d');
        case 'eCA_begin'
            par.CYeCA_b = sscanf(raw{i},'%*s %d');
        case 'eCA_end'
            par.CYeCA_e = sscanf(raw{i},'%*s %d');
        case 'iCA_begin'
            par.CYiCA_b = sscanf(raw{i},'%*s %d');
        case 'iCA_end'
            par.CYiCA_e = sscanf(raw{i},'%*s %d');
        case 'kuf'
            par.kuf = sscanf(raw{i},'%*s %e');
    end
end

%check to make sure all the fields were filled
fn = fieldnames(par);
for i = 1:length(fn)
    if isempty(par.(fn{i}))
        error('Error in load_params, missing parameter: %s\n',fn{i})
    end
end
        
%compute (or define) some additional parameters
par.TC = par.temp - 273.15;     %temp in celcius
par.sal = 35;                   %salinity assumed to be 35, doesn't have large effect on K1.                   
par.K1 = 10^-((3633.86./par.temp) - 61.2172 + 9.6777*log(par.temp) - 0.011555 .* par.sal + 0.0001152 .* (par.sal^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
par.K2 = 10^((-471.78/par.temp) - 25.929 + 3.16967 * log(par.temp) + 0.01781 * par.sal - 0.0001122*par.sal^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
par.Kw = exp(148.96502 - (13847.26 ./ par.temp) - 23.6521 .* log(par.temp) + (par.sal^.5).*((118.67 ./ par.temp) - 5.977 + 1.0495 .* log(par.temp)) - 0.01615 .* par.sal); %ion product of water, CO2 methods DOE 1994
par.O2sat = 230E-9;             %mol/cm3 at 20C in seawater
    
%diffusion coefficients
DCF = diff_coef(par.TC,par.sal, 1);        %diffusion coefficients (in cm2/s) at current temp, salinity, final variable is pressure: assume 1 atm.
par.Dc  = DCF(3);
par.Db  = DCF(8);       %diffusivity of HCO3- in cm2/s
    
%modify parameters for use in model
par.Vc = par.Ac * par.tx;                %volume of coral tissue
par.Vs  = (par.DBL .* 1) .* par.Ac;  %volume to consider as surface layer (5% of DBL thickness);
par.kur = par.kuf * (10^-par.pHe)./par.K1;              %uncatalyzed rate of HCO3 dehydration
par.bfrac = (1/ (1+ par.K2./(10^-par.pHe)));       % fraction of "h" pool that is HCO3-
    
return