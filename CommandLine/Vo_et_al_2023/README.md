# SSIT - Fitting and Experiment Design for HIV Reporter using MCP-GFP or smiFISH Labels

This folder contains codes needed to use the SSIT to fit models to experimental data and remake figures from the manuscript:

H. Vo et al, Analysis and design of single-cell experiments to harvest fluctuation information while rejecting measurement noise, 2023.

To run these codes, first make sure to follow the installation instructions for the general SSIT found in the README at: https://github.com/MunskyGroup/SSIT. In particular, you will need to install the Matlab Tensor Toolbox.

This folder contains two particular scripts of interest.  For each, it is recommended that you open these in the matlab editor and run one cell at a time.

"Vo_et_al_FittingExperimentalData.m" provides commands to run all data fitting routines. To complete fits, you will need to run the different sections multiple times, which will take a couple days of computing depending on your computer speed and number of available CPUs.

"Vo_et_al_PlottingResults.m" provides commands to make figures for the paper. 