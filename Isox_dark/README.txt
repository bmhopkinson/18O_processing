-This set of programs analyzes 18O exchange from CO2 catalyzed by cells in the dark. It calculates cellular transfer coefficients for CO2 (fc) and HCO3- (fb) , and the intracellular CO2 hydration rate (kcf), which is presumably catalyzed by carbonic anhydrase. 

-A sample data file is given: dark_test.txt, along with its parameter file: dark_test.par. To analyze this data call the isox_dark program at the Matlab command prompt:
>>isox_dark('dark_test',2)
provide the base name of the file without any file extensions. the option 2, describes the file format (see load_data.m for details).

- readouts from the optimization procedure should be shown as the program is running, converging after ~60 iterations. The output values should be kcf: 588 /s, fc: 2.7E-8 cm3/s, fb: 1.5E-11 cm3/s.  After completion of the parameter fitting, a plot of the raw data and model fit will be displayed. 

-a small output file is generated: dark_test_passive_params.out, which contains the fitted parameters. additionally a file called iso_exch_dark.out is written, which contains the model fit to the data. the parameters derived from this analysis are appended to the input parameter file (dark_test.par).


- the model used to fit the data is based on:
Tu C, Wynns GC, McMurray RE, and Silverman DN (1978) CO2 kinetics in red cell suspensions measured by 18O exchange. J. Biol. Chem. 253: 8178-8184.
with minor modifications and notation as discussed in:
Hopkinson BM, Dupont CL, Allen AE, and Morel FMM (2011) Efficiency of the CO2 concentrating mechanism of marine diatoms. PNAS. 

