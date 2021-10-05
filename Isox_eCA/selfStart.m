function kinit = selfStart(time,Obsdata,par)
%self starter function to determine good initial parameters for
%optimization
global FIG_NUM;
k = [1E-9 1E-8 5E-8 1E-7 5E-7 1E-6];      %fc values to test

imax = length(k);
RSS = zeros(1,imax);
par.return = 0;         %make sure eCA_lsq function returns only the predicted 13CO2 data
%determine residual sum of squares (RSS) for fits using each potential starting
%value

for i = 1:imax
    kvar = k(i);
    CO2fit = eCA_lsq(kvar,time,par);       %get predicted CO2 time course for kvar
    residuals = CO2fit - Obsdata;
    [m,n] = size(residuals);
    residuals = reshape(residuals',1,m*n);      %convert to vector
    RSS(1,i) = sum(residuals.^2);                %sum of squares
end

%find values with minimum RSS
[C,I] = min(RSS);

%optimization works better when starting at initial values below the actual
%values, so pick starting value one below that with minimium RSS
if I == 1
    kinit = k(I);
else
    kinit = k(I-1);
end


return
