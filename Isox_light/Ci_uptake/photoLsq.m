function residuals = photoLsq(x,tp,photo)

xint = round(x);       %convert to integer value
photos = zeros(size(photo));                        %initialize fitted photosynthesis vector
fit_photo1 = polyfit(tp(1:xint),photo(1:xint),1);       %fit line to ramping up portion
photos(1:xint) = polyval(fit_photo1,tp(1:xint));        %pull out values at particular timepoint
fit_photo2 = polyfit(tp(xint+1:end), photo(xint+1:end),1);  %fit line to steady state portion

%v1 = polyval(fit_photo2,tp(xint+1:end));    %pull out values at particular timepoints    
%v2 = photos(xint+1:end);
%size(v1)
%size(v2)

photos(xint+1:end) = polyval(fit_photo2,tp(xint+1:end));

residuals = photos - photo;     %modeled - observed

return;