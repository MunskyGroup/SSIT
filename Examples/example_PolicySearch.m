%% example_PolicySearch
% Simple policy search examples for control of two-state and three-state models. 
clear
addpath(genpath('../src'));

%% Example 1, a two-state model:  Transcription and Translation
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

%% Example 2, a three-state model:  DUSP1
% We use a simplified DUSP2 model where there is an upstream transcription
% factor (GR) that activates a gene.  Once active, the gene can transcribe
% nuclear RNA, which can later decay or leave the nucleus.
clear
DUSP1 = SSIT; 
DUSP1.species = {'offGene';'onGene';'rna'}; 
DUSP1.initialCondition = [2;0;0];         

% Define propensity functions and input signals:
DUSP1.propensityFunctions = {'kon*IGR*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
DUSP1.inputExpressions = {'IGR', ...
    '1 + a1*exp(-r1*(t-0)).*(t>=0)*(t<5) + a2*exp(-r2*(t-5)).*(t>=5)*(t<10) + a3*exp(-r3*(t-10)).*(t>=10)*(t<15) + a4*exp(-r4*(t-15)).*(t>=15)'}; 
DUSP1.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; 
DUSP1.parameters = ({'koff',0.014;'kon',0.002;'kr',1;'gr',0.004;...
                      'a1', 20; 'r1', 0.04; ...
                      'a2', 18; 'r2', 0.05; ...
                      'a3', 16; 'r3', 0.06; ...
                      'a4', 14; 'r4', 0.07}); 
DUSP1.tSpan = [0,5,10,15,20];
DUSP1.fspOptions.fspTol = 1e-5;
DUSP1 = DUSP1.formPropensitiesGeneral('DUSP1', true);
DUSP1.solutionScheme = 'FSP';
DUSP1.fspOptions.fspTol = 1e-7;
[~, DUSP1.fspOptions.bounds] = DUSP1.solve;
DUSP1.solutionScheme = 'fspSens';
DUSP1.fspOptions.fspTol = 1e-6;
OrigPars = [DUSP1.parameters{:,2}]';
tic
for i = 1:100
    DUSP1.parameters(:,2) = num2cell(OrigPars.*(1+0.1*randn(size(OrigPars))));
    SensSoln = DUSP1.solve;
    fspFIM = DUSP1.computeFIM(SensSoln.sens);
end
toc

%% DUSP1 Policy Search
DUSP1_PS = DUSP1;
DUSP1_PS.fspOptions.usePiecewiseFSP = true;
DUSP1_PS.solutionScheme = 'FSP';

logPolicy = zeros(1,4);
DUSP1_PS = DUSP1_PS.formPropensitiesGeneral('DUSP1_PS', true);

Target = [0,1,2,3,4];
loss = @(logPolicy)computeLossDUSP1(logPolicy,DUSP1_PS,Target);

tic
logPolicy = fminsearch(loss,logPolicy);
timeFindPolicy = toc
finalPolicy = exp(logPolicy)

[finalLoss,meanVal] = computeLossDUSP1(logPolicy,DUSP1_PS,Target);
finalLoss
plot(DUSP1_PS.tSpan,meanVal,DUSP1_PS.tSpan,Target,'x')


function [loss,meanVal] = computeLossDUSP1(logPolicy,DUSP1_PS,Target)
DUSP1_PS.parameters([5, 7, 9, 11],2) = num2cell(exp(logPolicy));
FSPsoln = DUSP1_PS.solve;
meanVal = zeros(1,5);
for i = 1:5   
    P = squeeze(sum(double(FSPsoln.fsp{i}.p.data),[1,2]));
    meanVal(i) = (0:length(P)-1)*P;
end
loss = sum((Target - meanVal).^2);
end
