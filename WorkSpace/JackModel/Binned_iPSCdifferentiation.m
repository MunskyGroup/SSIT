%% Using the SSIT to fit Multiple iPSC differentiation using machine learned groups
% scRNA seq of iPSC differentiation data is used to preform trajectory
% analysis using SSIT to get additional information to normal trajectory
% analysis that are unique to SSIT
close all 
clear
addpath(genpath('../src'));

%% Reading in the data
% goal is to make this universal for the groups because grouping may be
% changed later to betteraddpath('C:\Users\Jack\Documents\GitHub\SSIT_Base\src\CommandLine');addpath(genpath('../src')); fit biological groups instead of machine learned
% groups
dataloc = "Z:\home\formanj\scRNAseq_model\UpDownProject\Story\SSITData.csv";
data = readtable(dataloc);


%% Create Model
Model = SSIT;

% species
groupnames = data.Properties.VariableNames;
index = find(contains(groupnames, 'Group'));
speciesnames = cell(length(index),1);
for i = 1:length(speciesnames)
    speciesnames(i) = groupnames(index(i));
end

% speciesnames = {'Group_0'; 'Group_1'; 'Group_2'; 'Group_3'; 'Group_4'; 'Group_5'; ...
%  'Group_6'; 'Group_7'; 'Group_8'; 'Group_9';};

Model.species = speciesnames;

% initial conditions
ic = zeros([height(speciesnames), 1]);
ic(1) = 1;
Model.initialCondition = ic;

% Parameters, Propensities & Stoichiometery
pars = cell(height(speciesnames)*height(speciesnames), 2);
indx = 1;
propen = cell(height(speciesnames)*height(speciesnames), 1);
stoich = zeros([length(speciesnames), length(propen)]);
for i = 1:height(speciesnames)
    for j = 1:height(speciesnames)
        pars(indx, 1) = {'k_' + string(i-1) + string(j-1)};
        pars(indx, 2) = {1};
        propen(indx) = {string(speciesnames(i)) + '*' + pars(indx, 1)};
        stoich(i, indx) = -1;
        stoich(j, indx) = stoich(j, indx) + 1;
        indx = indx + 1;
    end
end
Model.parameters = pars;
Model.propensityFunctions = propen;
Model.stoichiometry = stoich;

Model.fspOptions.initApproxSS = true;
Model.dataSet = [];

%% Fit Model on all data
c = cell(length(speciesnames),2);
for i = 1:length(speciesnames)
    c{i, 1} = speciesnames(i);
    c{i, 2} = speciesnames(i);
end
Model = Model.loadData(dataloc, {'Group_0', 'Group_0'; ...
                                'Group_1', 'Group_1'; ...
                                'Group_2', 'Group_2'; ...
                                'Group_3', 'Group_3'; ...
                                'Group_4', 'Group_4'; ...
                                'Group_5', 'Group_5'; ...
                                'Group_6', 'Group_6'; ...
                                'Group_7', 'Group_7'; ...
                                'Group_8', 'Group_8'; ...
                                'Group_ 9', 'Group_9';});


ModelGR.fspOptions.initApproxSS = true;
ModelGR.dataSet = [];
%%
Model.solutionScheme = 'FSP';
[fspSoln,Model.fspOptions.bounds] = Model.solve;
Model.makePlot(fspSoln,'marginals')
Model.makeFitPlot

%%
Model.solutionScheme = 'FSP';
fitOptions = optimset('Display','iter','MaxIter',100);
fitOptions.suppressFSPExpansion = true; 
fitParameters = Model.maximizeLikelihood([],fitOptions);
Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(fitParameters);
Model.makeFitPlot

























