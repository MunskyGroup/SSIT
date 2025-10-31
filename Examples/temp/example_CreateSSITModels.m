%% example_CreateSSITModels
%  Example script to show how to create models in SSIT.
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
Model.parameters = ({'kon',30; 'koff',30; 'kr',100; 'gr',0.005});

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

% Update propensity function for the gene activation reaction:
STL1_Model.propensityFunctions{1} = 'offGene * IHog';

% Define the time-varying TF/MAPK input signal:
STL1_Model.inputExpressions = {'IHog',...
                               '(a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0))'};

% Add the new parameters from the TF/MAPK input signal
STL1_Model = STL1_Model.addParameter({'a0',0.01;'a1',1;'r1',0.4;'r2',.1});

% Print a summary of STL1_Model
STL1_Model.summarizeModel

%% Ex.(3) Load and modify a pre-existing SSIT model

% Several simple models are already coded in SSIT and can be quickly 
% loaded and modified. Below are currently available SSIT models (Mar2025). 

% ModelChoice = 'BirthDeath';      % One species problem
% ModelChoice = 'CentralDogma':    % Two species problem
% ModelChoice = 'BurstingGene':    % Two species problem
% ModelChoice = 'ToggleSwitch';    % Two species problem 
%                                    (non-linear toggle switch)
% ModelChoice = 'ToggleSwitch2';   % Two species problem
% ModelChoice = 'CentralDogmaTV';  % Two species problem (mRNa and protein) 
%                                    with time-varying transcription rate
% ModelChoice = 'Repressilator';               % Three species problem
% ModelChoice = 'BurstingSpatialCentralDogma'; % Four species problem 
%                                            (gene state, nuclear mRNa, 
%                                            cytoplasmic mRNA, and protein)
ModelChoice = 'RepressilatorGenes'; % Nine species problem

% Let's load one of the pre-existing SSIT models
RepGenes_Model = SSIT(ModelChoice);

% View a summary of the 'RepressilatorGenes' model
RepGenes_Model.summarizeModel

               