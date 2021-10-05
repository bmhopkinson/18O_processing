-This set of programs analyzes the uncatalyzed rate of 18O exchange from CO2 prior to the addition of cells or carbonic anhydrase. It calculates the uncatalyzed rate of CO2 hydration in the solution (kuf).

-A sample data file is given: uncat_test.txt, along with its parameter file: uncat_test.par. To analyze this data call the isox_uncat program at the Matlab command prompt:
>>isox_uncat('uncat_test'),2
provide the base name of the file without any file extensions.the option 2, describes the file format (see load_data.m for details).

- readouts from the optimization procedure should be shown as the program is running, converging after ~5 iterations. After completion of the parameter fitting, a plot of the raw data and model fit will be displayed. 

- the fitted kuf is appended to the input parameter file for using in further processes (such as isox_dark.m).

- the equations to fit the data in uncat_lsq.m are from:
Silverman DN (1974) A new approach to measuring the rate of rapid bicarbonate exchange across membranes. Mol. Pharmacol. 10: 820-836. 