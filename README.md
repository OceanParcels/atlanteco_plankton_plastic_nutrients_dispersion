# Source code for article titled: Separation timescales of vertically migrating zooplankton and other (a)biotic materials in the Benguela system

We use GLOB16 high resolution ocean model data to study separation timescales between different particle types in the Benguela system. The code used for this study is organized as following: 

```data/```: Files generated for release locations of particles in the Benguela system and global set for daily dawn and dusk times.

###### Programs
Code for the analysis can be used in the following order-

1. ```src/processing/extract_Benguela_region.ipynb```: To create the set of release locations across the Benguela system.
2. ```src/processing/SunriseSunsetTool.py```: Compute dawn and dusk times for a global grid of 2°x2°for each day of the year.
3. ```src/simulations/Benguela_fwd_simulation.py```: Run simulations for particles with diffferent properties: 2D, 3D, DVM and sinking.
4. ```src/analysis/FullDistanceDaysAnalysis.py```: Code to compute the separation distance between particle pairs and CDF with separation distance of 100 km.

###### Figures
Plots in the manuscript have been generated using the following codes:
- Figure 1: ```src/visualizations/Domain_plot.ipynb```: Domain plot
- Figure 3, 4 & 5: ```src/analysis/Separation_Distance_analysis.ipynb```
- Figure 6, 7B, 8 & 9: ```src/analysis/Mesoscale_separation_timescales.ipynb```
- Figure A1: ```src/processing/DVM_kernel_test.ipynb```
- Figure A2, A3 & A4:```src/analysis/Separation_Distance_analysis.ipynb```
- Figure A5: ```src/analysis/Sensitivity_analysis.ipynb```
