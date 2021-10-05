function Ciup = load_Ciup(infile)
%loads photosynthesis and net Ci uptake (CO2 and HCO3-) data.

Ciup = struct('tc',[],'O2',[],'totCO2',[],'Photo',[],'Cup',[],'Bup',[]);

fid = fopen(infile,'r');

line = fgetl(fid);      %discard header line

i = 1;
while ~feof(fid)
    line = fgetl(fid);
    A = sscanf(line, '%f %e %e %e %e %e');
    tc(1,i) = A(1);
    O2(1,i) = A(2);
    totCO2(1,i) = A(3);
    Photo(1,i) = A(4);
    Cup(1,i) = A(5);
    Bup(1,i) = A(6);
    i=i+1;
end

Ciup.tc = tc;
Ciup.O2 = O2;
Ciup.totCO2 = totCO2;
Ciup.Photo = Photo;
Ciup.Photo(1,end)
Ciup.Cup = Cup;
Ciup.Bup = Bup;

return