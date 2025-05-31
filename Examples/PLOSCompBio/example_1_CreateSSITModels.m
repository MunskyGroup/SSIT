%% example_1_CreateSSITModels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.1: Creating, Saving, and Loading Models in SSIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries:
clear
close all
addpath(genpath('../../src'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Create a simple Bursting Gene model in SSIT
%%  This model consists of 3 species: 
%   an inactive gene ('offGene'), an activated gene ('onGene'), and mRNA. 
%%  There are four reactions, each with a unique rate parameter: 
%   activation ('kon'), inactivation ('koff'), transcription ('kr'), 
%   and mRNA degradation ('gr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create an SSIT instance and call it 'Model':
Model = SSIT;    

% Set species names for bursting gene expression model:
Model.species = {'offGene';'onGene';'mRNA'}; 

% Set initial condition:
Model.initialCondition = [2;0;0];           

% Set stoichiometry of reactions:
Model.stoichiometry = [-1,1,0,0;...
                        1,-1,0,0;...
                        0,0,1,-1]; 

% Set propensity functions:
Model.propensityFunctions = {'kon * offGene';'koff * onGene';...
                             'kr * onGene';'gr * mRNA'}; 

% Set initial guesses for parameters:
Model.parameters = ({'kon',0.2; 'koff',0.2; 'kr',10; 'gr',5});

% Print a summary of STL1 Model:
Model.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Copy the Bursting Gene Model above to a new model, 'STL1'.  
% The species, reactions, and stochiometries remain the same; 
% however, a change is made to the propensity for the first reaction 
% (gene activation) so that its rate is controlled by a time-varying 
% MAPK input signal. As a result, new parameters must also be added.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a copy of the simple Bursting Gene model from above:
STL1 = Model;

% Update propensity function for the gene activation reaction:
STL1.propensityFunctions{1} = 'offGene * IHog';

% Define the time-varying TF/MAPK input signal:
STL1.inputExpressions = {'IHog',...
                              '(a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0))'};

% Add the new parameters from the TF/MAPK input signal:
STL1.parameters = ({'koff',0.2; 'kr',10; 'gr',5;...
                    'a0',5; 'a1',10; 'r1',0.004; 'r2',.01});

% Print a summary of STL1 Model:
STL1.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex.(3) Load and modify a pre-set SSIT model

% Several simple models are already coded in SSIT and can be quickly 
% loaded and modified. Below are currently available SSIT models (2025). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Let's load one of the pre-existing SSIT models:
RepGenes_Model = SSIT(ModelChoice);

% View a summary of the 'RepressilatorGenes' model:
RepGenes_Model.summarizeModel