%% Using the SSIT to fit Multiple iPSC differentiation using machine learned groups
% scRNA seq of iPSC differentiation data is used to preform trajectory
% analysis using SSIT to get additional information to normal trajectory
% analysis that are unique to SSIT
close all 
clear
addpath(genpath('../../src'));

%% Reading in the data
% goal is to make this universal for the groups because grouping may be
% changed later to betteraddpath('C:\Users\Jack\Documents\GitHub\SSIT_Base\src\CommandLine');addpath(genpath('../src')); fit biological groups instead of machine learned
% groups
dataloc = "SSITData.csv";
data = readtable(dataloc);


%% Create Model
Model = SSIT;

% species
varNames = data.Properties.VariableNames;
% index = find(contains(varNames, 'leiden'));
% unique(data.leiden)
% speciesnames = cell(unique(data.leiden),1);
% for i = 1:length(speciesnames)
%     speciesnames(i) = groupnames(index(i));
% end

% speciesnames = {'Group_0'; 'Group_1'; 'Group_2'; 'Group_3'; 'Group_4'; 'Group_5'; ...
%  'Group_6'; 'Group_7'; 'Group_8'; 'Group_9';};

Model.species = {'Group'};

% nGroups = max(data.leiden);
nGroups = max(data.leiden);

Pars = rand(max(data.leiden)+1); Pars = Pars-diag(diag(Pars));

% initial conditions
% ic = zeros([height(speciesnames), 1]);
% ic(1) = 1;
Model.initialCondition = 1;

pars = {};
n =0;
for i = 0:nGroups
    for j = 0:nGroups
        if i~=j
            n = n+1;
            pars = [pars;{['k',num2str(i),'_',num2str(j)],Pars(i+1,j+1)}];
            % pars = [pars;{['k',num2str(n)],Pars(i+1,j+1)}];
        end
    end
end
   
reactions = cell(nGroups+nGroups,1);
stoichs = zeros(1,nGroups+nGroups);
stoichs(1:nGroups) = [1:nGroups];
stoichs(nGroups+1:end) = -[1:nGroups];

for i = 0:nGroups-1
    for j = i+1:nGroups
         reactions{j-i} =         [reactions{j-i},        '+k',num2str(i),'_',num2str(j),'*(Group==',num2str(i),')'];
         reactions{nGroups+j-i} = [reactions{nGroups+j-i},'+k',num2str(j),'_',num2str(i),'*(Group==',num2str(j),')'];
         % reactions{j-i} =         [reactions{j-i},        '+k',num2str(n),'*(Group==',num2str(i),')'];
         % reactions{nGroups+j-i} = [reactions{nGroups+j-i},'+k',num2str(n),'*(Group==',num2str(j),')'];
    end
end

Model.parameters = pars;
Model.propensityFunctions = reactions;
Model.stoichiometry = stoichs;
Model.formPropensitiesGeneral('JackPropens')

return
%%

% cell(height(speciesnames)*height(speciesnames), 2);
% indx = 1;
% propen = cell(height(speciesnames)*height(speciesnames), 1);
% stoich = zeros([length(speciesnames), length(propen)]);
% for i = 1:height(speciesnames)
%     for j = 1:height(speciesnames)
%         pars(indx, 1) = {'k_' + string(i-1) + string(j-1)};
%         pars(indx, 2) = {1};
%         propen(indx) = {string(speciesnames(i)) + '*' + pars(indx, 1)};
%         stoich(i, indx) = -1;
%         stoich(j, indx) = stoich(j, indx) + 1;
%         indx = indx + 1;
%     end
% end
% Model.parameters = pars;
% Model.propensityFunctions = propen;
% Model.stoichiometry = stoich;
% 
% Model.fspOptions.initApproxSS = true;
% Model.dataSet = [];

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

























