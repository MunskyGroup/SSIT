
% clear all
close all
parfor i=4:6
    switch i
        case 1
            iterativeExperimentRunner('Poisson','simulated','FIMopt',10,1,i)
        case 2
            iterativeExperimentRunner('Poisson','simulated','FIMopt',10,1,i)
        case 3
            iterativeExperimentRunner('Poisson','simulated','uniform',10,1,i)
        case 4            
            iterativeExperimentRunner('DUSP1','simulated','FIMopt',10,1,i)
        case 5
            iterativeExperimentRunner('DUSP1','simulated','random',10,1,i)
        case 6
            iterativeExperimentRunner('DUSP1','simulated','uniform',10,1,i)
    end
end

%% Poisson Results
clear det*
load IterativeExperimentResults_Poisson_simulated_FIMopt_2
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_FIM(i-1) = det(covLogMH{i});
    detFIMInv_FIM(i-1) = predictCov(FIMcurrentExptSaved{i});
    detFIMTrueInv_FIM(i-1) = predictCov(FIMcurrentExptTrueSaved{i});
end

load IterativeExperimentResults_Poisson_simulated_uniform_2
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Unif(i-1) = det(covLogMH{i});
    detFIMInv_Unif(i-1) = predictCov(FIMcurrentExptSaved{i});
    detFIMTrueInv_Unif(i-1) = predictCov(FIMcurrentExptTrueSaved{i});
end

load IterativeExperimentResults_Poisson_simulated_random_2
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Rand(i-1) = det(covLogMH{i});
    detFIMInv_Rand(i-1) = predictCov(FIMcurrentExptSaved{i});
    detFIMTrueInv_Rand(i-1) = predictCov(FIMcurrentExptTrueSaved{i});
end

figure(1); clf;
plot(2:nExpt,detCov_FIM,'b',...
    2:nExpt,detCov_Unif,'r',...
    2:nExpt,detCov_Rand,'m',...
    'linewidth',2)

hold on
plot(2:nExpt,detFIMInv_FIM,'--b',...
    2:nExpt,detFIMInv_Unif,'--r',...
    2:nExpt,detFIMInv_Rand,'--m',...
    'linewidth',2)

plot(2:nExpt,detFIMTrueInv_FIM,'-.b',...
    2:nExpt,detFIMTrueInv_Unif,'-.r',...
    2:nExpt,detFIMTrueInv_Rand,'-.m',...
    'linewidth',2)

set(gca,"FontSize",16,'yscale','log')
legend('FIM','Uniform','Random')

%% DUSP1 Results
clear det*
load IterativeExperimentResults_DUSP1_simulated_FIMopt_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_FIM(i-1) = det(covLogMH{i});
    detFIMInv_FIM(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]);
    detFIMTrueInv_FIM(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]);
end

load IterativeExperimentResults_DUSP1_simulated_uniform_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Unif(i-1) = det(covLogMH{i});
    detFIMInv_Unif(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]);
    detFIMTrueInv_Unif(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]);
end

load IterativeExperimentResults_DUSP1_simulated_random_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Rand(i-1) = det(covLogMH{i});
    detFIMInv_Rand(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]);
    detFIMTrueInv_Rand(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]);
end

figure(2); clf;
plot(2:nExpt,detCov_FIM,'b',...
    2:nExpt,detCov_Unif,'r',...
    2:nExpt,detCov_Rand,'m',...
    'linewidth',2)

hold on
plot(2:nExpt,detFIMInv_FIM,'--b',...
    2:nExpt,detFIMInv_Unif,'--r',...
    2:nExpt,detFIMInv_Rand,'--m',...
    'linewidth',2)

plot(2:nExpt,detFIMTrueInv_FIM,'-.b',...
    2:nExpt,detFIMTrueInv_Unif,'-.r',...
    2:nExpt,detFIMTrueInv_Rand,'-.m',...
    'linewidth',2)
    
set(gca,"FontSize",16,'yscale','log')
% legend('MH','Predicted','Exact')
% legend('FIM','Uniform','Random')

function predictedDetCov = predictCov(fimSet,parSet)
arguments
    fimSet
    parSet = [1:length(fimSet{1})];
end
detCovPreds = zeros(1,length(fimSet));
for i = 1:length(fimSet)
    detCovPreds(i) = det(fimSet{i}(parSet,parSet)^(-1));
end
predictedDetCov = mean(detCovPreds);
end