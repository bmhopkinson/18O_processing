function [time, O2, Ar, CO2, C13O2] = load_data(data_file)
%load MIMS data from text file
%file format:  time, H2O, O2, Ar, CO2, 13CO2, 13C18O16O,
%13C28O2 (by increasing mass)

%read in MS data from the rest of the file    
fdat = fopen(data_file,'r');                      %skip header line of file

  C = textscan(fdat, '%u %f %f %f %f %f %f %f %f');   %read data into a cell array
  time  = C{2}';                      %time in seconds
  O2    = C{4}';                      %mass 36, Oxygen
  Ar    = C{5}';                      %mass 40, Argon
  CO2   = C{6}';                      %mass 44, CO2
  C13O2 = [C{7}';C{8}';C{9}'];        %masses 45, 47, 49 13CO2s  

fclose(fdat);
end