-This program processes 18O exchange data collected on corals in the absence and presence of a eCA inhibitor (e.g. DBAZ) to calculate surface eCA activity. It also fits internal CA activity and membrane permeabilities to data obtained the presence of the eCA inhibitor, but these parameters are not well constrained and should only be used to account for the effects of iCA activity earlier in the run (prior to addition of DBAZ).

- a sample data file is give: eCA_test.txt, along with its parameter file eCA_test.par. To analyze this data call the isox_coral_eCA program at the Matlab command prompt.
>> isox_coral_eCA('eCA_test')
- the program will attempt to optimize the model fit to the data and should converge after ~20 iterations, giving approximate values of: ksf = 1.7, kcf = 15, Pc = 1.1E-2, Pb = 3.3E-4.
- you may need to modify the initial conditions for the optimization to get good fits (lowest f(x) in the iterations readout, which is the error function being optimized). To do so open "fit_coral_eCA.m" and modify the "kvar_init" values. 
- This program requires the Matlab Optimization Toolbox and the Statistics Toolbox