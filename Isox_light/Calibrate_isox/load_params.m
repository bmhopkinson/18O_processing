function [par, init] = load_params(file)
%opens parameter file and loads in data 

fid = fopen(file,'r');

ct = 1;
raw={};
while ~feof(fid)
    raw{ct}=fgetl(fid);
    ct=ct+1;
end
fclose(fid);

par = struct('cycles',[],'temp',[],'K1',[],'evol',[],'cells',[],'cyvol',[],'xvol',[],'cySA',[],'xSA',[],'kuf',[],'kif',[],'Pc',[],'Pb',[],...
    'pHe',[],'DIC',[],'CYback_b',[],'CYback_e',[],'CYuncat_b',[],'CYuncat_e',[],'CYdark_b',[],'CYdark_e',[],'CYlight_b',[],'CYlight_e',[],...
    'kur',[],'kir',[]);
init = struct('be',[],'ci',[],'bi',[]);

for i = 1:length(raw)
    [id] = sscanf(raw{i},'%s[\t]%*s');
    
    switch id
        
        case 'cycles'
            par.cycles =sscanf(raw{i},'%*s %d');
        case 'DIC'
            par.DIC = sscanf(raw{i},'%*s %f');
        case 'temp'
            par.temp = 273.15 + sscanf(raw{i},'%*s %f');
        case 'evol'
            par.evol = sscanf(raw{i},'%*s %f');
        case 'enrich_factor'
            par.enrich = sscanf(raw{i},'%*s %f');
        case 'cells/mL'
            par.cells = sscanf(raw{i},'%*s %e');
        case 'cell_vol'
            par.cyvol = sscanf(raw{i},'%*s %e');
        case 'chl_vol'
            par.xvol =  sscanf(raw{i},'%*s %e');
        case 'cell_SA'
            par.cySA = sscanf(raw{i},'%*s %e');
        case 'chl_SA'
            par.xSA = sscanf(raw{i},'%*s %e');
        case 'kuf'
            par.kuf = sscanf(raw{i},'%*s %e');
        case 'kcf'
            par.kif = sscanf(raw{i},'%*s %e');
        case 'Pc'
            par.Pc = sscanf(raw{i},'%*s %e');
        case 'Pc_x'
            par.Pc_x = sscanf(raw{i},'%*s %e');
        case 'Pb'
            par.Pb =  sscanf(raw{i},'%*s %e');
        case 'pHe'
            par.pHe = sscanf(raw{i},'%*s %f');
        case 'HCO3_ext'
            init.be = sscanf(raw{i},'%*s %e %e %e %e');
        case 'CO2_in'
            init.ci = sscanf(raw{i},'%*s %e %e %e');
        case 'HCO3_in'
            init.bi = sscanf(raw{i},'%*s %e %e %e %e');
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
    end
end
    par.sal = 35;                   %salinity assumed to be 35, doesn't have large effect on K1.                   
    par.K1 = 10^-((3633.86./par.temp) - 61.2172 + 9.6777*log(par.temp) - 0.011555 .* par.sal + 0.0001152 .* (par.sal^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
    par.K2 = 10^((-471.78/par.temp) - 25.929 + 3.16967 * log(par.temp) + 0.01781 * par.sal - 0.0001122*par.sal^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
    par.Kw = exp(148.96502 - (13847.26 ./ par.temp) - 23.6521 .* log(par.temp) + (par.sal^.5).*((118.67 ./ par.temp) - 5.977 + 1.0495 .* log(par.temp)) - 0.01615 .* par.sal); %ion product of water, CO2 methods DOE 1994
    
    
    %modify parameters for use in model
    par.xvol = (par.cells * par.xvol)/10;                        %total volume of chloroplasts in MIMS chamber
    par.xSA = (par.cells * par.xSA);
    par.cyvol = (par.cells * par.cyvol) - par.xvol;         %total cytosolic volume in MIMS chamber
    par.cySA = par.cells * par.cySA;
    fb =(1/ (1+ par.K2./(10^-par.pHe)));                    % fraction of "h" pool that is HCO3- to correct HCO3 dehydration rates.                                                        % corrected for fraction of "HCO3-" in model which is actul
    par.kur = par.kuf * fb * (10^-par.pHe)./par.K1;              % uncatalyzed rate of HCO3 dehydration,
    par.kir = par.kif * fb * (10^-par.pHe)./par.K1;              % catalyzed rate of HCO3 dehydration
    
    par;
    init;
    
    return