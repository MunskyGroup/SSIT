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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Copy the 'STL1' model for 4-state dynamics.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a copy of the STL1 model from above:
STL1_4state = STL1;

% Set species names for STL1_4state:
STL1_4state.species = {'s1';'s2';'s3';'s4';'mRNA'};

% Set initial condition:
STL1_4state.initialCondition = [2;0;0;0;0];  

% Keep the time-varying TF/MAPK input signal:
STL1_4state.inputExpressions = {'IHog',...
                              '(a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0))'};

% Set stoichiometry of reactions:
STL1_4state.stoichiometry = [0,-1,1,0,0,0,0,0,0,0,0;...  % gene state 1
                             0,1,-1,0,-1,1,0,0,0,0,0;... % gene state 2
                             0,0,0,0,1,-1,0,-1,1,0,0;... % gene state 3
                             0,0,0,0,0,0,0,1,-1,0,0;...  % gene state 4
                             1,0,0,1,0,0,1,0,0,1,-1];    % mRNA
                % Reactions: 1,2,3,4,5,6,7,8,9,10,11

% Set propensity functions:
STL1_4state.propensityFunctions = {'kr1*s1';'k12*s1';'k21*s2*IHog';...
                                   'kr2*s2';'k23*s2';'k32*s3';...
                                   'kr3*s3';'k34*s3';'k43*s4';...
                                   'kr4*s4';'dr * mRNA'}; 

% Add the new parameters for the 4 state model:
STL1_4state.parameters = ({'k12',0.2; 'k23',0.2; 'k34',0.2;...
                          'k21',0.2; 'k32',0.2; 'k43',0.2; ...
                          'kr1',10; 'kr2',10; 'kr3',10; 'kr4',10; ...
                          'dr',5; 'a0',5; 'a1',10; 'r1',0.004; 'r2',.01});

% Print a summary of STL1 Model:
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex.(4) Load and modify a pre-set SSIT model

% Several simple models are already coded in SSIT and can be quickly 
% loaded and modified. Below are currently available SSIT models (2025). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ModelChoice = 'BirthDeath';      % One species problem
% ModelChoice = 'CentralDogma';    % Two species problem
% ModelChoice = 'BurstingGene';    % Two species problem
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