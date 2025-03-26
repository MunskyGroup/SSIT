%% example_CreateSSITModels
%  Example script to show how to create a model in SSIT.
clear
close all
addpath(genpath('../src'));

%% Ex.(1): Create a simple Bursting Gene model in SSIT

%% This model consists of 3 species: 
% an inactive gene ('offGene'), an activated gene ('onGene'), and mRNA. 

%% There are four reactions, each with a unique rate parameter: 
% activation ('kon'), inactivation ('koff'), transcription ('kr'), 
% and mRNA degradation ('gr')

% Create an SSIT instance and call it 'Model':
Model = SSIT;    

% Set species names for bursting gene expression model:
Model.species = {'offGene';'onGene';'mRNA'}; 

% Set stoichiometry of reactions:
Model.stoichiometry = [-1,1,0,0;...
                        1,-1,0,0;...
                        0,0,1,-1]; 

% Set propensity functions:
Model.propensityFunctions = {'kon * offGene';'koff * onGene';...
                             'kr * onGene';'gr * mRNA'}; 

% Set initial guesses for parameters:
Model.parameters = ({'kon',30; 'koff',100; 'kr',0.005; 'gr',0.01});

% Set initial condition (one 'offGene'):
Model.initialCondition = [1;0;0]; 

% Print a summary of our Model:
Model.summarizeModel


%% Ex.(2) Create an SSIT model for real, time-varying experimental data: 
%%           yeast STL1 data provided by the Vanderbilt Q_BIO Group

% Create, solve, and fit a CME model to single-cell smFISH data. 
% The dataset is from yeast STL1 collected in Dr. Gregor Neuert's 
% laboratory at Vanderbilt University.

% Copy the Bursting Gene Model above to a new model, 'STL1_Model'.  
% The species, reactions, and stochiometries remain the same; 
% however, a change is made to the propensity for the first reaction 
% (gene activation) so that its rate is controlled by a time-varying 
% MAPK input signal. As a result, new parameters must also be added.

% Create a copy of the simple Bursting Gene model from above:
STL1_Model = Model;

% Update propensity functions:
STL1_Model.propensityFunctions = {'offGene * IHog'; 'koff * onGene';...
                                    'kr * onGene'; 'gr * mRNA'};

% Define the time-varying TF/MAPK input signal:
STL1_Model.inputExpressions = {'IHog',...
                               '(a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0))'};

% Add the new parameters from the TF/MAPK input signal
STL1_Model = STL1_Model.addParameter({'a0',0.01;'a1',1;'r1',0.4;'r2',.1});

% Print a summary of STL1_Model
STL1_Model.summarizeModel