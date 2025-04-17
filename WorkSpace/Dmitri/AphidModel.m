close all
clear all
addpath(genpath('../../src'));

Model = SSIT;
% Model.species = {'aphid';'predator'};
Model.species = {'aphid'};
Model.initialCondition = [0];

% Model.propensityFunctions = {'k1','g1*aphid','k2','g2*predator'};
Model.propensityFunctions = {'k1+k2*aphid^eta/(M^eta+aphid^eta)','g1*aphid'};

% Model.stoichiometry = [1,-1,0,0;...
%     0,0,1,-1];
Model.stoichiometry = [1,-1];

% Model.parameters = ({'k1',5;'g1',0.5;'k2',2;'g2',0.3});
Model.parameters = ({'k1',5;'g1',0.5;'k2',10;'eta',2;'M',100});

Model.fspOptions.initApproxSS = false;
Model.fittingOptions.modelVarsToFit = (1:2);

Model = Model.formPropensitiesGeneral('DmitriAphids');

[FSPGrSoln,Model.fspOptions.bounds] = Model.solve;

%%
% Model = Model.loadData("AphidsNewFormatB.csv",...
%     {'aphid','sca'});
Model = Model.loadData("AphidsNewFormatB.csv",...
    {'aphid','sca'},...
    {'plant_date','E';...
    'seed_trtmt','U';...
    'spray_trtmt','USP'});

% {'Condition','GR_timesweep';'Dex_Conc',GRfitCases{i,2}}
%%
Model.tSpan = unique(Model.dataSet.times);
Model.initialTime = Model.tSpan(1);

Model.makeFitPlot
set(gca,'xlim',[Model.tSpan(1),Model.tSpan(end)])

%%
% Model.parameters(1:2,2) = {2;0.1};
Model.fspOptions.bounds(2) = 1000;

Model.initialCondition = [0];
fitOptions = optimset('Display','iter','MaxIter',300);
pars = [Model.parameters{:,2}];
for i = 1:1
    fitOptions.suppressFSPExpansion = true;
    pars = Model.maximizeLikelihood(...
        pars, fitOptions);
    Model.parameters(:,2) = num2cell(pars);
end
Model.makeFitPlot
