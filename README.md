# MATLAB Stochastic System Identification Toolkit (SSIT) for modeling single-cell fluorescence microscopy data

![SSIT](https://github.com/MunskyGroup/SSIT/blob/main/images/GraphicalAbstract.png)

Authors: Huy Vo, Joshua Cook, Brian Munsky

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The SSIT allows users to specifiy and solve the chemical master equation for discreate stochastic models, especially those used for the analysis of single-cell gene regulaton.  

To learn more about the FSP theory that underlies the SSIT, please see the slides from our [Nov. 3, 2022 BPPB Seminar](https://github.com/MunskyGroup/SSIT/blob/main/images/BPPBSeminarNov2022.pdf)


The SSIT includes command line tools and a graphical user interface to:
- Build, save, and load models
-- Create from propensity functions and stoichiometries
-- Load from / export to SimBiology
-- Load from / export to SBML
- Solve models
-- Solve using ODE analyses and basic moment closure analyses
-- Generate synthetic data from models using Stochastic Simulations (parallel)
-- Solve CME directly using the Finite State Projection algorithm, included automated FSP state set selection/expansion
-- All methods support non-linear, time varying propensity functions, including logical statements 
- Compute sensitivity of solutions to parameter variations
- Compute first passage or escape time distributions for complex trajectories
- Load experimental single-cell data (e.g., from processes smFISH images)
- Fit Models to Experimental data
-- Compute the likelihood of data given model
-- Maximize likelihood using gradient and non-gradient based searches
-- Run efficient Metropolis Hastings algorithm (with custom proposal distribution or proposal distribution based on FIM) to estimate parameter uncertainties given single-cell data
-- Include custom priors on parameter distributions for Bayesian analysis
- Model and reject measurement noise
-- Calibrate empirical probability distortion operator (PDO) to quantify effects of image distortion
-- Include image distortions in parameter estimation
- Improve Experiment Designs
-- Compute the Fisher Information Matrix for CME models 
-- Sample over model uncertainty (e.g., model prior or posterior from previous experiment round) for iterative experiment design.
-- Search experiment design space to find optimally informative experiments
-- Automatically adjust designs to account for image distortion effects
- Explore effects of Extrinsic Noise in parameters
- Form reduce order models
-- Run hybrid models with deterministic and stochastic species
-- Reduce models using Quasi-Steady Approximations on Fast Species
-- Reduce models using Eigenvalue decomposition
-- Reduce models using coarse meshes
-- Reduce models using Principle Orthogonal Decomposition
- Compare multiple models to different data sets with shared parameter sets
-- Identify parameters that change with genetic/environmental/experimental conditions
- Many Examples

# Dependencies
For all basic functionalities:
- MATLAB R2021b or later.
- Symbolic Computing Toolbox.
- Global Optimization Toolbox (for model fitting only)
- Parallel Computing Toolbox (optional). 
- SimBiology Toolbox (for loading/saving SBML models only)

# Installation
Clone this package to a local folder on your computer. Then add the path to that folder (with subfolders) into MATLAB's search path. You can then call all functions from MATLAB. 

# Testing
To test your installation, navigate to the folder SSIT/tests and run the following test routines:
PoissonTest % Tests various solution schemes, data generation/loading and model parameterization for a 1-species model
Poisson2DTest % Tests various solutions for a 2-species model
PoissonTVTest % Tests various solutions for a 1-species Time-varying model
multiModelTests % Tests various solutions for a combinations of multiple models with different data sets.
modelReductionTests % Tests various model reduction schemes for more efficient solutions of the Chemical Master Equation
miscelaneousTests % Tests other aspects, including GUI functionality, loading/saving to SBME and SimBiology

# Acknowledgements

The SSIT tools in this repository make use of sparse tensors using the Tensor Toolbox for MATLAB (version 3.2.1) provided by Brett W. Bader, Tamara G. Kolda and others at www.tensortoolbox.org under   Users of this software should cite the TTM creators at:

* B. W. Bader and T. G. Kolda, Efficient MATLAB Computations with Sparse and Factored Tensors, SIAM J. Scientific Computing, 30(1):205-231, 2007, http://dx.doi.org/10.1137/060676489. 

The provided SSIT tools also make use a modified version of Expokit for the solution of time-invariant master equations (although these codes are not used for the current publication). Users of this software should cite the creators at:

* Sidje, R. B., Expokit Software Package for Computing Matrix Exponentials, ACM Trans. Math. Softw., 24:1, 1998.

# Getting Started
The SSIT provides two basic interaction options: (1) command line tools and (2) a graphical user interface.

## GUI Version  
A GUI version of the SSIT has much of the functionality, and is a great way to familiarize yourself with the approach.  However, for intenssive research tasks, we strongly recommend using the command line tools. To get started with the GUI, compile and launch to tool kit with the following commands:

>> src2app;

>> A = SSIT_GUI;

You should then see the model loading and building page of the graphical interface, and you are off to the races...
![SSIT](https://github.com/MunskyGroup/SSIT/blob/main/images/SSITGUI_snapshot.png)

## Command Line Version

To get started with the Command line Tools, navigate to the directory "CommandLine" and open one of the tutorial scripts "example_XXX.m".  Or you caan start creating and solving models as follows.

Example for generating an FSI model and fitting it to smFISH data for Dusp1 activation following glucocorticoid stimulation:

Define SSIT Model

>> Model = SSIT;

>> Model.species = {'OnGene';'rna'};

>> Model.initialCondition = [0; 0];

>> Model.propensityFunctions = {'kon * IGR * (2-OnGene)'; 'koff * OnGene'; 'kr * OnGene'; 'gr * rna'};

>> Model.stoichiometry = [1,-1,0,0; 0,0,1,-1];

>> Model.inputExpressions = {'IGR','a0 + a1 * exp(-r1 * t) * (1-exp(-r2 * t)) * (t>0)'};

>> Model.parameters = ({'koff',0.14; 'kon',0.14; 'kr',25; 'gr',0.01; 'a0',0.006; 'a1',0.4; 'r1',0.04; 'r2',0.1});
    
>> Model.fspOptions.initApproxSS = true;  % Model is assumed to start at steady state at t=0;

Load and Fit smFISH Data
    
>> Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'rna','RNA_nuc'});
    
>> Model.tSpan = unique([Model.initialTime,Model.dataSet.times]);
>> fitOptions = optimset('Display','iter','MaxIter',100);
>> pars,likelihood] = Model.maximizeLikelihood([],fitOptions);

Update Model and Make Plots of Results

>> Model.parameters(:,2) = num2cell(pars);

>> Model.makeFitPlot

You should arrive at a fit of the model to the experimentally measured Dusp1 mRNA distributions looking something like this:

![SSIT](https://github.com/MunskyGroup/SSIT/blob/main/images/Dusp1Fit.png)






