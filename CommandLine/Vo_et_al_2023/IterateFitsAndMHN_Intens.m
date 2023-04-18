% Iterate on cluster
function IterateFitsAndMHN_Intens(iJob)

vars.mhScaling = 0.8;

if iJob<=3
    fileStr = 'FISHTrue_2StateBurst_4Pars_Prior2x'
    kJob = iJob
    vars.priorScale=2;
    vars.modelChoice = '2stateBurst';
    vars.modelVarsToFit = [2:5];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(7.12),log10(0.012)]';
    vars.sigLog10Prior = [1,1,.5,.5,.5]'*vars.priorScale;
%     vars.mhScaling = [0.8,0.8,0.4,0.4,0.4];
    vars.initialFspBounds = [0 0 0 0 3 3 3 1200];
elseif iJob<=7
    fileStr = 'FISHTrue_2StateBurst_4Pars_Prior1x'
    kJob = iJob-4
    vars.modelChoice = '2stateBurst';
    vars.priorScale=1;
    vars.modelVarsToFit = [2:5];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(7.12),log10(0.012)]';
    vars.sigLog10Prior = [1,1,.5,.5,.5]'*vars.priorScale;
%     vars.mhScaling = [0.8,0.8,0.4,0.4,0.4];
    vars.initialFspBounds = [0 0 0 0 3 3 3 1200];
elseif iJob==8
    fileStr = 'FISHTrue_3Pars_TwoStatePoisson'
    vars.modelChoice = '2statePoisson';
    kJob = 3
    vars.priorScale=2;
    vars.modelVarsToFit = [2:4];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(0.012)]';
    vars.sigLog10Prior = [1,1,3,.5]'*vars.priorScale;
    vars.initialFspBounds = [0 0 0 3 3 1200];
elseif iJob==18
    fileStr = 'FISHTrue_3Pars_TwoStatePoisson'
    vars.modelChoice = '2statePoisson';
    kJob = 2
    vars.priorScale=2;
    vars.modelVarsToFit = [2:4];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(0.012)]';
    vars.sigLog10Prior = [1,1,3,.5]'*vars.priorScale;
    vars.initialFspBounds = [0 0 0 3 3 1200];
elseif iJob==19
    fileStr = 'FISHTrue_3Pars_TwoStatePoisson'
    vars.modelChoice = '2statePoisson';
    kJob = 1
    vars.priorScale=2;
    vars.modelVarsToFit = [2:4];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(0.012)]';
    vars.sigLog10Prior = [1,1,3,.5]'*vars.priorScale;
    vars.initialFspBounds = [0 0 0 3 3 1200];
elseif iJob==102
    fileStr = 'FISHTrue_3Pars_TwoStatePoisson_Guess'
    vars.modelChoice = '2statePoisson';
    kJob = 3
    vars.priorScale=2;
    vars.modelVarsToFit = [1:4];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(0.012)]';
    vars.sigLog10Prior = [1,1,3,.5]'*vars.priorScale;
    vars.initialFspBounds = [0 0 0 3 3 1200];
elseif iJob==9
    fileStr = 'FISHTrue_5Pars_TwoStateBursts'
    vars.modelChoice = '2stateBurst';
    kJob = 3
    vars.priorScale=2;
    vars.modelVarsToFit = [1:5];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(7.12),log10(0.012)]';
    vars.sigLog10Prior = [1,1,.5,.5,.5]'*vars.priorScale;
    vars.initialFspBounds = [0 0 0 0 3 3 3 1200];
elseif iJob==10
    fileStr = 'FISHTrueThreeState_6pars'
    vars.modelChoice = '3state';
    kJob = 3
    vars.priorScale=2;
    vars.modelVarsToFit = [1:6];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(7.12),log10(0.012),3]';
    vars.sigLog10Prior = [1,1,.5,.5,.5,1]'*vars.priorScale;
    vars.initialFspBounds = [0 0 0 0 3 3 3 1200];
elseif iJob==101
    fileStr = 'FISHTrueThreeState_6pars_Guess'
    vars.modelChoice = '3state';
    kJob = 3
    vars.priorScale=2;
    vars.modelVarsToFit = [1:6];
    vars.muLog10Prior = [-4,-4,log10(0.2),log10(7.12),log10(0.012),3]';
    vars.sigLog10Prior = [1,1,.5,.5,.5,1]'*vars.priorScale;
    vars.initialFspBounds = [0 0 0 0 3 3 3 1200];
elseif iJob==20
    fileStr = 'FISHTrue_OneStateBurst'
    kJob = 3
    vars.priorScale=1;
    vars.modelChoice = '1stateBurst';
    vars.modelVarsToFit = [1:3];
    vars.muLog10Prior = [log10(0.2),log10(7.12),log10(0.012)]';
    vars.sigLog10Prior = [.5,.5,.5]'*vars.priorScale;
    vars.mhScaling = 0.2;
end

if iJob<=7
    vars.otherParGuesses = {'FISHTrue_2StateBurst_4Pars_Prior1x_0_18',...
        'FISHTrue_2StateBurst_4Pars_Prior1x_0_300',...
        'FISHTrue_2StateBurst_4Pars_Prior1x_0_18_300',...
        'FISHTrue_2StateBurst_4Pars_Prior2x_0_18',...
        'FISHTrue_2StateBurst_4Pars_Prior2x_0_300',...
        'FISHTrue_2StateBurst_4Pars_Prior2x_0_18_300'};
elseif iJob==8||iJob==18||iJob==19
    vars.otherParGuesses = {'FISHTrue_3Pars_TwoStatePoisson_0_300',...
        'FISHTrue_3Pars_TwoStatePoisson_0_18_300',...
        'FISHTrue_3Pars_TwoStatePoisson_0_18'};
end

try
    parpool(5);
    parfevalOnAll(@warning,0,'off','all')
catch
    parfevalOnAll(@warning,0,'off','all')
end

TMP = SSIT;
clear TMP;
try
    addpath(genpath('../../MatlabToolboxes/tensor_toolbox-v3.2.1/'))
catch
end

vars.doFit = 1;
vars.display = 'none';
vars.iter = 10;
vars.fitCases = [1:5];

switch kJob
    case 0
        vars.timeSet = 0; %Data to fit: {'0','0_18','0_300','0_18_300'}
    case 1
        vars.timeSet = [0,18]; %Data to fit: {'0','0_18','0_300','0_18_300'}
        fileStr=[fileStr,'_0_18'];
    case 2
        vars.timeSet = [0,300]; %Data to fit: {'0','0_18','0_300','0_18_300'}
        fileStr=[fileStr,'_0_300'];
    case 3
        vars.timeSet = [0,18,300]; %Data to fit: {'0','0_18','0_300','0_18_300'}
        fileStr=[fileStr,'_0_18_300'];
end
    
FittingFunctionsFISHTrue(1,fileStr,vars)

for i=1:5
    if i==5
        vars.nMH=20000;
        vars.nFIMsamples = 20;
    else
        vars.nMH=1000;
        vars.nFIMsamples = 3;
    end

    for iStep = [32,3,4,5]
        iStep
        FittingFunctionsFISHTrue(iStep,fileStr,vars)
    end

end

i = 1;
while exist([fileStr,'_',num2str(i),'.mat'],'file')
    i=i+1;
end