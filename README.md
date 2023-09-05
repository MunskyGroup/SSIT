# MATLAB Stochastic System Identification Toolkit (SSIT) for modeling single-cell fluorescence microscopy data

![SSIT](https://github.com/MunskyGroup/SSIT/blob/main/images/GraphicalAbstract.png)

Authors: Huy Vo, Joshua Cook, Brian Munsky

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The SSIT allows users to specifiy and solve the chemical master equation for discreate stochastic models, especially those used for the analysis of single-cell gene regulaton.  

To learn more about the FSP theory that underlies the SSIT, please see the slide from our [Nov. 3, 2022 BPPB Seminar](https://github.com/MunskyGroup/SSIT/blob/main/images/BPPBSeminarNov2022.pdf)


The SSIT includes command line tools and a graphical user interface to:
- Build, save, and load models
- Generate synthetic data from models using Stochastic Simulations
- Solve models using the Finite State Projection algorithm
- Compute sensitivity of FSP solutions to parameter variations
- Load experimental smFISH data
- Compute/Maximize the likelihood of data givien model
- Run Metropolis Hastings algorithm to estimate parametr uncertainties given single-cell data
- Compute the Fisher Information Matrix for CME models
- Search experiment design space to find optimally informative experiments

# Dependencies
For all basic functionalities:
- MATLAB R2021b or later.
- Symbolic Computing Toolbox.
- Global Optimization Toolbox.
- Parallel Computing Toolbox. 
- [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org/).  You will need to make sure to add the TTB to the Matlab path before running the SSIT.

# Installation
Clone this package to a local folder on your computer. Then add the path to that folder (with subfolders) into MATLAB's search path. You can then call all functions from MATLAB. 

# Acknowledgements

The SSIT tools in this repository make use of sparse tensors using the Tensor Toolbox for MATLAB (version 3.2.1) provided by Brett W. Bader, Tamara G. Kolda and others at www.tensortoolbox.org under   Users of this software should cite the TTM creators at:

* B. W. Bader and T. G. Kolda, Efficient MATLAB Computations with Sparse and Factored Tensors, SIAM J. Scientific Computing, 30(1):205-231, 2007, http://dx.doi.org/10.1137/060676489. 

The provided SSIT tools also make use a modified version of Expokit for the solution of time-invariant master equations (although these codes are not used for the current publication). Users of this software should cite the creators at:

* Sidje, R. B., Expokit Software Package for Computing Matrix Exponentials, ACM Trans. Math. Softw., 24:1, 1998.

# Getting Started
The SSIT provides two basic interaction options: (1) command line tools and (2) a graphical user interface.

## GUI Version  
To get started with the GUI, compile and launch to tool kit with the following commands:

>> src2app;

>> A = SSIT_GUI;

You should then see the model loading and building page of the graphical interface, and you are off to the races...
![SSIT](https://github.com/MunskyGroup/SSIT/blob/main/images/SSITGUI_snapshot.png)

## Command Line Version

To get started with the Command line Tools, navigate to the directory "CommandLine" and open one of the tutorial scripts "example_XXX.m".  Or you caan start creating and solving models as follows.

Example for generating an FSI model and fitting it to smFISH data for Dusp1 activation following glucocorticoid stimulation:

Define SSIT Model

>> Model = SSIT;

>> Model.species = {'x1';'x2'};

>> Model.initialCondition = [0; 0];

>> Model.propensityFunctions = {'kon * IGR * (2-x1)'; 'koff * x1'; 'kr * x1'; 'gr * x2'};

>> Model.stoichiometry = [1,-1,0,0; 0,0,1,-1];

>> Model.inputExpressions = {'IGR','a0 + a1 * exp(-r1 * t) * (1-exp(-r2 * t)) * (t>0)'};

>> Model.parameters = ({'koff',0.14; 'kon',0.14; 'kr',25; 'gr',0.01; 'a0',0.006; 'a1',0.4; 'r1',0.04; 'r2',0.1});
    
>> Model.initialTime = -120;  % large negative time to simulate steady state at t=0

Load and Fit smFISH Data
    
>> Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});
    
>> Model.tSpan = unique([Model.initialTime,Model.dataSet.times]);
>> fitOptions = optimset('Display','iter','MaxIter',100);
>> pars,likelihood] = Model.maximizeLikelihood([],fitOptions);

Update Model and Make Plots of Results

>> Model.parameters(:,2) = num2cell(pars);

>> Model.makeFitPlot

You should arrive at a fit of the model to the experimentally measured Dusp1 mRNA distributions looking something like this:

![SSIT](https://github.com/MunskyGroup/SSIT/blob/main/images/Dusp1Fit.png)






