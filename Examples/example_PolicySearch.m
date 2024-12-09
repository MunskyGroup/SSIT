%% example_PolicySearch
% Simple policy search example for control of two-state model. 
clear
addpath(genpath('../src'));

%% A two-state model:  Transcription and Translation
% First create the SSIT model 
TT = SSIT;
TT.species = {'rna','protein'};
TT.initialCondition = [1;0];
TT.propensityFunctions = {'kr'; 'gr*rna'; 'kp*rna'; 'gp*protein'};
TT.stoichiometry = [1,-1,0,0; 0,0,1,-1];
TT.parameters = {'kr',0.5; 'gr',0.1; 'kp',0.5; 'gp',0.1}; 
TT.tSpan = [0,5,10,15];
TT.fspOptions.fspTol = 1e-5;
TT = TT.formPropensitiesGeneral('TT', true);
TT.solutionScheme = 'FSP';
TT.fspOptions.fspTol = 1e-7;
[~, TT.fspOptions.bounds] = TT.solve;
TT.solutionScheme = 'fspSens';
TT.fspOptions.fspTol = 1e-6;

OrigPars = [TT.parameters{:,2}]';

for i = 1:100
    TT.parameters(:,2) = num2cell(OrigPars.*(1+0.1*randn(size(OrigPars))));
    SensSoln = TT.solve;
    fspFIM = TT.computeFIM(SensSoln.sens);
end


%% Transcription and Translation Policy Search
TT_PS = TT;
TT_PS.fspOptions.usePiecewiseFSP = true;
TT_PS.solutionScheme = 'FSP';
TT_PS.parameters = {'kr',0.5; 'gr1',0.1; 'gr2',0.1; 'gr3',0.1; 'kp',0.5; 'gp',0.1}; 
TT_PS.propensityFunctions = {'kr'; 'gr1*rna*(t<5)+gr2*rna*(t>=5)*(t<10)+gr3*rna*(t>=10)'; 'kp*rna'; 'gp*protein'};
TT_PS = TT_PS.formPropensitiesGeneral('TT', true);

Target = [0,1,2,3]; 

logPolicy = zeros(1,3);
loss = @(logPolicy)computeLoss(logPolicy,TT_PS,Target);

logPolicy = fminsearch(loss,logPolicy);
finalPolicy = exp(logPolicy)

[finalLoss,meanVal] = computeLoss(logPolicy,TT_PS,Target);
plot(TT_PS.tSpan,meanVal,TT_PS.tSpan,Target,'x')


function [loss,meanVal] = computeLoss(logPolicy,TT_PS,Target)
TT_PS.parameters(2:4,2) = num2cell(exp(logPolicy));
FSPsoln = TT_PS.solve;
    meanVal = zeros(1,4); 
    for i = 1:4
        P = sum(double(FSPsoln.fsp{i}.p.data),1);
        meanVal(i) = (0:length(P)-1)*P(:);
    end
    loss = sum((Target - meanVal).^2);
end
