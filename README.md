# MATLAB Stochastic System Identification Toolkit (SSIT) for modeling single-cell fluorescence microscopy data

## Dependencies
For all basic functionalities:
- MATLAB R2021b or later.
- Symbolic Computing Toolbox.
- Global Optimization Toolbox.
- Parallel Computing Toolbox.
- [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org/).

Optionally, Sundial's CVODE and CVODES can be used to solve the reduced CME (optional). This requires that [Sundials](https://computing.llnl.gov/projects/sundials) library is installed and is in the computer's search path.

## Installation
Clone this package to a local folder on your computer. Then add the path to that folder (with subfolders) into MATLAB's search path. You can then call all functions from MATLAB. 

## Getting Started
The SSIT provides two basic interaction options: (1) command line tools and (2) a graphical user interface.

To get started with the Command line Tools, navigate to the directory "CommandLine" and open one of the tutorial scripts "example_XXX.m".

To get started with the GUI, compile and launch to tool kit with the following commands:

%: src2app;

%: A = SSIT_GUI;
