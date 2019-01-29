function p = load_params(file, index)
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
p = struct('DIC',[],'Ve',[],'pHe',[],'temp',[],'enrich',[],'cells',[],'Vcell',[],'kuf',[],'kcf',[],'fc',[],'fb',[],...
    'Ac',[],'tx',[],'Vhom',[],'Vadd',[],'ksf',[],...
    'CYback_b',[],'CYback_e',[],'CYuncat_b',[],'CYuncat_e',[],'CYhom_b',[],'CYhom_e',[],'CYinh_b',[],'CYinh_e',[]);

for i = 1:length(raw)
    [id] = sscanf(raw{i},'%s[\t]%*s');
    
    switch id
        
        case 'DIC'
            p.DIC = (sscanf(raw{i},'%*s %f'))*1E-9;       %DIC in mol/cm3
        case 'evol'
            p.Ve = sscanf(raw{i},'%*s %f');       %assay volume cm3
        case 'pHe'
            p.pHe = sscanf(raw{i},'%*s %f');      %pH of assay solution
        case 'temp'
            p.temp = 273.15 + sscanf(raw{i},'%*s %f');        %temp in K
        case 'enrich_factor'
            p.enrich = sscanf(raw{i},'%*s %f');               %accounts for enrichment of Ci prior to dilution by sample addition
        case 'residual_zoox'
            p.cells = sscanf(raw{i},'%*s %e');    %concentration of residual zoox in assay solution cells/cm3
        case 'cell_vol'
            p.Vcell = sscanf(raw{i},'%*s %e');    %volume of a single zoox cell
        case 'kuf'
            p.kuf = sscanf(raw{i},'%*s %e');      %uncatalyzed CO2 hydration rate /s
        case 'kcf'
            p.kcf = sscanf(raw{i},'%*s %e');      %internal CA activity of zoox /s
        case 'fc'
            p.fc = sscanf(raw{i},'%*s %e');       %cellular transfer coeff for CO2 cm3/s
        case 'fb'
            p.fb =  sscanf(raw{i},'%*s %e');      %cellular transfer coeff for HCO3- cm3/s
        case 'SA'
            p.Ac = sscanf(raw{i},'%*s %e');       %surface area of coral removed cm2
        case 'tissue_thickness'
            p.tx = sscanf(raw{i},'%*s %e');       %coral effective tissue thickness cm
        case 'homg_vol'
            p.Vhom = sscanf(raw{i},'%*s %e');     %homogenate volume cm3
        case 'addn_vol'
            p.Vadd = sscanf(raw{i},'%*s %e');     %volume of sample homogenate added to assay cm3
        case 'ksf'
            p.ksf = sscanf(raw{i},'%*s %e');      %eCA activity of coral cm/s
        case 'background_begin'
            p.CYback_b = sscanf(raw{i},'%*s %d');
            p.CYback_b = p.CYback_b - index(1) + 1;     %correct if cycle numbers don't start at 1 
        case 'background_end'
            p.CYback_e = sscanf(raw{i},'%*s %d');
            p.CYback_e = p.CYback_e - index(1) + 1;
        case 'uncat_begin'
            p.CYuncat_b = sscanf(raw{i},'%*s %d');
            p.CYuncat_b = p.CYuncat_b - index(1) + 1;
        case 'uncat_end'
            p.CYuncat_e = sscanf(raw{i},'%*s %d');
            p.CYuncat_e = p.CYuncat_e - index(1) + 1;
        case 'homog_begin'
            p.CYhom_b = sscanf(raw{i},'%*s %d');
            p.CYhom_b = p.CYhom_b - index(1) + 1;
        case 'homog_end'
            p.CYhom_e = sscanf(raw{i},'%*s %d');
            p.CYhom_e = p.CYhom_e - index(1) + 1;
        case 'inhib_begin'
            p.CYinh_b = sscanf(raw{i},'%*s %d');
            p.CYinh_b = p.CYinh_b - index(1) + 1;
        case 'inhib_end'
            p.CYinh_e = sscanf(raw{i},'%*s %d');
            p.CYinh_e = p.CYinh_e - index(1) + 1;

    end
end

%check to make sure all the fields were filled
fn = fieldnames(p);
for i = 1:length(fn)
    if isempty(p.(fn{i}))
        error('Error in load_params, missing parameter: %s\n',fn{i})
    end
end
        
%compute (or define) some additional parameters
    p.sal = 35;                   %salinity assumed to be 35, doesn't have large effect on K1.                   
    p.K1 = 10^-((3633.86./p.temp) - 61.2172 + 9.6777*log(p.temp) - 0.011555 .* p.sal + 0.0001152 .* (p.sal^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
    p.K2 = 10^((-471.78/p.temp) - 25.929 + 3.16967 * log(p.temp) + 0.01781 * p.sal - 0.0001122*p.sal^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
    p.Kw = exp(148.96502 - (13847.26 ./ p.temp) - 23.6521 .* log(p.temp) + (p.sal^.5).*((118.67 ./ p.temp) - 5.977 + 1.0495 .* log(p.temp)) - 0.01615 .* p.sal); %ion product of water, CO2 methods DOE 1994
    p.O2sat = 230E-9;             %mol/cm3 at 20C in seawater
    p.h = 10^-p.pHe;            %H+ conc
    
    %modify parameters for use in model
    p.Vc = p.Ac * p.tx;                %volume of coral tissue
    p.bfrac = (1/ (1+ p.K2./(10^-p.pHe)));       % fraction of "h" pool that is HCO3-
    p.kur = p.bfrac .* p.kuf .* (p.h)./p.K1;              %uncatalyzed rate of HCO3 dehydration
    p.kcr = p.bfrac .* p.kcf .* (p.h)./p.K1;    %catalyzed rate of HCO3- dehydration in zoox cells
    p.ksf_s = (p.ksf .* p.Ac .* p.Vadd ./ (p.Vhom .* p.Ve));        %coral eCA activity when homogenized (/s)
    p.ksr_s = p.bfrac .* p.ksf_s .* p.h./p.K1; 
   
    
    p;
return