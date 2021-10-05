function CO2_sum = apply_cfs(cfs, x, params)

cfs(1,1:3) = cfs;

CO2_sigs = params;

CALIB = diag(cfs);                  % create a matrix w/ calibration factors on the diagonal 

CO2_calib = CALIB * CO2_sigs;       % apply input calibration factors to CO2 current data

CO2_sum = sum(CO2_calib,1);         

return