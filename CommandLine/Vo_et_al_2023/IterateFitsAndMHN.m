% Iterate on cluster
function IterateFitsAndMHN(iJob)
if iJob<=3
    fileStr = 'Fit_using_SS_Priorx_FixedDeg'
    kJob = iJob
    vars.priorScale=1;
    vars.modelVarsToFit = [1:4];
elseif iJob<=7
    fileStr = 'Fit_using_SS_Prior1x'
    kJob = iJob-4
    vars.priorScale=1;
    vars.modelVarsToFit = [1:5];
elseif iJob<=11
    fileStr = 'Fit_using_SS_Prior4x'
    kJob = iJob-8
    vars.priorScale=4;
    vars.modelVarsToFit = [1:5];
elseif iJob<=15
    fileStr = 'FixedOnRate_using_SS_Prior1x'
    kJob = iJob-12
    vars.priorScale=1;
    vars.modelVarsToFit = [2:5];
elseif iJob<=19
    fileStr = 'FixedOnRate_using_SS_Prior4x'
    kJob = iJob-16
    vars.priorScale=4;
    vars.modelVarsToFit = [2:5];
elseif iJob<=23
    fileStr = 'FixedOnRate_PDOon0_300_Prior1x'
    kJob = iJob-20
    vars.priorScale=1;
    vars.modelVarsToFit = [2:5];
    vars.pdoTimes = [0,300];
end

if iJob==0
    vars.mhScaling = 0.8;
elseif iJob==4
    vars.mhScaling = 0.6;
else
    vars.mhScaling = 0.8;
end

if iJob==21
    vars.mhScaling = 1;
end

try
    parpool(5);
catch
end
TMP = SSIT;
clear TMP;
try
    addpath(genpath('../../MatlabToolboxes/tensor_toolbox-v3.2.1/'))
catch
end

vars.doFit = 1;
vars.display = 'final';
vars.iter = 400;

switch kJob
    case 0
        vars.timeSet = '0';
    case 1
        vars.timeSet = '0_18';
        fileStr=[fileStr,'_0_18'];
    case 2
        vars.timeSet = '0_300';
        fileStr=[fileStr,'_0_300'];
    case 3
        vars.timeSet = '0_18_300';
        fileStr=[fileStr,'_0_18_300'];
end

Q = load(fileStr,'FIMZero','FIMZeroLog',...
    'ModelZero','chainResults','covLog','fimResults','lTrue2FISH',...
    'lTrue2MCP','fimResultsOne','parsets');
save([fileStr,'_BAK'],'Q')

FittingFunctionsCoLocalized(3,fileStr,vars)
FittingFunctionsCoLocalized(4,fileStr,vars)


for i=1:5
    if i==5
        vars.nMH=10000;
        vars.nFIMsamples = 20;
    else
        vars.nMH=1000;
        vars.nFIMsamples = 5;
    end
    switch kJob
        case 0
            FittingFunctionsCoLocalized(2,fileStr,vars)
        case 1
            FittingFunctionsCoLocalized(12,fileStr,vars)
        case 2
            FittingFunctionsCoLocalized(22,fileStr,vars)
        case 3
            FittingFunctionsCoLocalized(32,fileStr,vars)
    end
    FittingFunctionsCoLocalized(3,fileStr,vars)
    FittingFunctionsCoLocalized(4,fileStr,vars)
    FittingFunctionsCoLocalized(5,fileStr,vars)
end

i = 1;
while exist([fileStr,'_',num2str(i),'.mat'],'file')
    i=i+1;
end
load(fileStr,'FIMZero','FIMZeroLog',...
    'ModelZero','chainResults','covLog','fimResults','lTrue2FISH',...
    'lTrue2MCP','fimResultsOne','parsets');
save([fileStr,'_',num2str(i)],'FIMZero','FIMZeroLog',...
    'ModelZero','chainResults','covLog','fimResults','lTrue2FISH',...
    'lTrue2MCP','fimResultsOne','parsets')
