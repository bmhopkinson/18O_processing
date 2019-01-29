This program calculates internal carbonic anhydrase activity (iCA) by analysis of 18O exchange rates as catalyzed by coral tissue homogenate.
To determine iCA accurately, the program needs to know the coral eCA activity (it can also account for zooxanthellae iCA contamination, but this should not be an issue if zooxanthellae are removed by centrifugation).
The model is described in Hopkinson et al. 2015 J. Exp Biol "Internal carbonic anhydrase activity in the tissue of scleractinian corals is sufficient to support proposed roles in photosynthesis and calcification".

- a sample data file is give: coral_iCA_test.txt, along with its parameter file coral_iCA_test.par. To analyze this data call the isox_coral_iCA program at the Matlab command prompt.
>> isox_coral_iCA('coral_iCA_test')
- the program will attempt to optimize the model fit to the data and should converge after ~6 iterations, giving approximate values of: kif = 121 +- 3 /s
- you may need to modify the initial conditions for the optimization to get good fits (lowest f(x) in the iterations readout, which is the error function being optimized). To do so open "fit_coral_iCA.m" and modify the "kvar_init" value. 
- This program requires the Matlab Optimization Toolbox and the Statistics Toolbox