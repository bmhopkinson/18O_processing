
### <sup>18</sup>O_processing: Matlab scripts to analyze <sup>18</sup>O-CO<sub>2</sub> exchange data
This respository contains Matlab scripts to infer carbonic andhydrase activity, cell membrane permeability, and inorganic carbon fluxes from  <sup>18</sup>O-isotope exchange experiments.
Each folder contains scripts specialized for different conditions or organisms. Consult readmes and referenced papers in the individual folders for details.
In general, the experiments are conducted by adding <sup>18</sup>O-enriched inorganic carbon to a vessel containing phytoplankton, corals, or other organisms. 
The exchange of <sup>18</sup>O for <sup>16</sup>O from water is then monitored by Membrane Inlet Mass Spectrometry (MIMS). Carbonic anhydrase within the organisms accelerates <sup>18</sup>O removal.
The observed temporal dynamics of <sup>18</sup>O removal are then fit using a system of differntial equations to extract information on carbonic andhydrase activity, cell membrane permeability, and inorganic carbon fluxes.

The figure below shows an example from work on phytoplankton. The diagram on the left shows the inorganic carbon forms and fluxes under consideration in the model. 
The plot on the right shows observed <sup>18</sup>O removal dynamics (open circles) and the model fit (solid lines) from a marine diatom. 
![18O_Example](doc/18O_model_fit_example.jpg)

The associated system of differential equations is:
![18O_DEs](doc/18O_model_DEs.jpg)

see Hopkinson et al. 2011 PNAS 108: 3830-3837 for details.
