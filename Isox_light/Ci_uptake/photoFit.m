function photos = photoFit(tp, photo);
%model photosynthesis data as two linear regions: a ramp up and then
%approximate steady state
%optimize to determine break between these two regions

x0 = 20;        %inital guess for break between the ramp up and steady state
lb = 5;         %lower bound
ub = 50;        %upper bound

f=@(x)photoLsq(x,tp,photo);     %anonymous functio to pass time and data vector to fitting function

%set options for lsqnonlin
options = optimset('DiffMinChange',1,'Display','final','TolFun',1E-21);


[x, resnorm,residual,exitflag] = lsqnonlin(f,x0,lb,ub,options);

%apply optimzed x to determine fitted photosynthetic rates
xint = round(x);        %convert to integer value
photos = zeros(size(photo));                        %initialize fitted photosynthesis vector
fit_photo1 = polyfit(tp(1:xint),photo(1:xint),1);       %fit line to ramping up portion
photos(1:xint) = polyval(fit_photo1,tp(1:xint));        %pull out values at particular timepoint
fit_photo2 = polyfit(tp(xint+1:end), photo(xint+1:end),1);  %fit line to steady state portion
photos(xint+1:end) = polyval(fit_photo2,tp(xint+1:end));    %pull out values at particular timepoints    

photos;

return
