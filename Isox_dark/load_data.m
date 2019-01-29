function [time, O2, Ar, CO2, C13O2] = load_data(data_file, FILETYPE)
%load MIMS data,
%file format 1 = line number, time, H2O, CO2, 13CO2, 13C18O16O, 13C28O2, Ar, O2
%file format 2 = line number, time, H2O, O2, Ar, CO2, 13CO2, 13C18O16O,
%13C28O2 (by increasing mass)

%read in MS data from the rest of the file    
fdat = fopen(data_file,'r');                      %skip header line of file

if FILETYPE == 1				%old file mass orders
  C = textscan(fdat, '%u %f %f %f %f %f %f %f %f %f');  %read data into a cell array
  time  = C{2}';                      %time in seconds
  CO2   = C{4}';                      %mass 44 C12O16O16
  C13O2 = [C{5}';C{6}';C{7}'];        %mass 45, 47, 49 13CO2s
  Ar    = C{8}';                      %mass 40, Argon
  O2    = C{9}';                      %mass 32, Oxygen
    
elseif FILETYPE == 2			%new file mass orders, in order of increasing mass
  C = textscan(fdat, '%u %f %f %f %f %f %f %f %f');   %read data into a cell array
  time  = C{2}';                      %time in seconds
  O2    = C{4}';                      %mass 36, Oxygen
  Ar    = C{5}';                      %mass 40, Argon
  CO2   = C{6}';                      %mass 44, CO2
  C13O2 = [C{7}';C{8}';C{9}'];        %masses 45, 47, 49 13CO2s  
end
    
fclose(fdat);
end