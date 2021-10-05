### Overview
This set of programs analyzes 18O exchange from CO2 catalyzed by cells in the abscence of an eCA inhibitor to determine extracellular carbonic anhydrase activity.
The CO2 and HCO3- mass transfer coefficients and iCA activity must already be known and input in the .par file. These parameters can be determined by analysis of 18O-exchange catalyzed by cells in the prescence of an eCA inhibitor (isox_dark.m). 

### Usage
A sample data file is given: eCA_test.txt, along with its parameter file: eCA_test.par. To analyze this data call the isox_eCA program at the Matlab command prompt:
  ```Shell
  isox_eCA('eCA_test')
  ```

provide the base name of the files without any file extensions.

the program will access parameters from the .par file and data from the .txt file, which must have the same base name. For information on the file formats see the sample files, templates in isox_dark code, and the respective loading subroutines (load_params.m, load_data.m).
readouts from the optimization procedure should be shown as the program is running, converging after ~12 iterations. The fitted eCA value should be ~1.95E-6.  After completion of the parameter fitting, a plot of the raw data and model fit will be displayed. The fitted eCA activity is appened to the .par file.

two output files are generated: basename_eCA_fit.out, which contains the model fit, and basename_eCA_CO2.out, which contains the calibrated CO2 data. additionally a file called iso_exch_dark.out is written, which contains the model fit to the data.


### Citation
the model used to fit the data is from:
Hopkinson BM, Meile C, and Shen C (2013) Quantification of extracellular carbonic anhydrase activity in two marine diatoms and investigation of its role Plant Physiology 162: 1142-1152  doi: http:/​/​dx.​doi.​org/​10.​1104/​pp.​113.​217737 
