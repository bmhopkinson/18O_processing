function par = load_params(file)
%opens parameter file and loads in data 

fid = fopen(file,'r');

ct = 1;
raw={};
while ~feof(fid)
    raw{ct}=fgetl(fid);
    ct=ct+1;
end
fclose(fid);

%define all required parameters to be read from file in a structure and set as empty. later will
%check if any are left empty 
par = struct('DIC',[],'evol',[],'pHe',[],'temp',[],'cells',[],'cyvol',[],'cySA',[],'enrich',[],...
    'CYback_b',[],'CYback_e',[],'CYuncat_b',[],'CYuncat_e',[],'CYdark_b',[],'CYdark_e',[],'CYlight_b',[],...
    'CYlight_e',[]);

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

    end
end

%check to make sure all the fields were filled
    fn = fieldnames(par);
    for i = 1:length(fn)
        if isempty(par.(fn{i}))
            error('Error in load_params, missing parameter: %s\n',fn{i})
        end
    end

    %add some additional parameters
    par.sal = 35;                   %salinity assumed to be 35, doesn't have large effect on K1.                   
    par.K1 = 10^-((3633.86./par.temp) - 61.2172 + 9.6777*log(par.temp) - 0.011555 .* par.sal + 0.0001152 .* (par.sal^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
    par.K2 = 10^((-471.78/par.temp) - 25.929 + 3.16967 * log(par.temp) + 0.01781 * par.sal - 0.0001122*par.sal^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
    par.Kw = exp(148.96502 - (13847.26 ./ par.temp) - 23.6521 .* log(par.temp) + (par.sal^.5).*((118.67 ./ par.temp) - 5.977 + 1.0495 .* log(par.temp)) - 0.01615 .* par.sal); %ion product of water, CO2 methods DOE 1994
        
    par;
    
    return