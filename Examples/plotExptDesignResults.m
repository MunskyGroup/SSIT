%% Section to run all Experiment Designs (FIMopt, Random, Uniform, and Intuition)
% Each Design will save the results in a mat file that will be loaded into
% the next sections and plotted.
clear all
close all

nRound = 10; % Number of rounds for each experiment design
parfor i=4:7
    switch i
        case 1
            iterativeExperimentRunner('Poisson','simulated','FIMopt',nRound,1,i)
        case 2
            iterativeExperimentRunner('Poisson','simulated','FIMopt',nRound,1,i)
        case 3
            iterativeExperimentRunner('Poisson','simulated','uniform',nRound,1,i)
        case 4            
            iterativeExperimentRunner('DUSP1','simulated','FIMopt',nRound,1,i)
        case 5
            iterativeExperimentRunner('DUSP1','simulated','random',nRound,1,i)
        case 6
            iterativeExperimentRunner('DUSP1','simulated','uniform',nRound,1,i)
         case 7
            iterativeExperimentRunner('DUSP1','simulated','intuition',nRound,1,i)
    end
end


%% Load and Plot Poisson Results
clear det*

load IterativeExperimentResults_Poisson_simulated_FIMopt_2
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_FIM(i-1) = det(covLogMH{i}); % |COV| of MH samples
    detFIMInv_FIM(i-1) = predictCov(FIMcurrentExptSaved{i}); % Predicted |FIM^-1| 
    detFIMTrueInv_FIM(i-1) = predictCov(FIMcurrentExptTrueSaved{i}); % True |FIM^-1| 
end

load IterativeExperimentResults_Poisson_simulated_uniform_2
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Unif(i-1) = det(covLogMH{i}); % |COV| of MH samples
    detFIMInv_Unif(i-1) = predictCov(FIMcurrentExptSaved{i}); % Predicted |FIM^-1| 
    detFIMTrueInv_Unif(i-1) = predictCov(FIMcurrentExptTrueSaved{i}); % True |FIM^-1|
end

load IterativeExperimentResults_Poisson_simulated_random_2
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Rand(i-1) = det(covLogMH{i}); % |COV| of MH samples
    detFIMInv_Rand(i-1) = predictCov(FIMcurrentExptSaved{i}); % Predicted |FIM^-1| 
    detFIMTrueInv_Rand(i-1) = predictCov(FIMcurrentExptTrueSaved{i}); % True |FIM^-1|
end

figure(1); clf;
% Invisible plot lines for last legend entries
plot(2:nExpt,detCov_FIM,'k',...
    2:nExpt,detFIMInv_FIM,'--k',...
    2:nExpt,detFIMTrueInv_FIM,'-.k',...
    'linewidth',2);

% PLot |COV| of each Design
hold on
plot(2:nExpt,detCov_FIM,'b',....
    2:nExpt,detCov_Unif,'r',...
    2:nExpt,detCov_Rand,'m',...
    'linewidth',2);

% Plot predicted FIM of each Design
plot(2:nExpt,detFIMInv_FIM,'--b',...
    2:nExpt,detFIMInv_Unif,'--r',...
    2:nExpt,detFIMInv_Rand,'--m',...
    'linewidth',2);

% Plot rue |FIM| of each Design 
plot(2:nExpt,detFIMTrueInv_FIM,'-.b',...
    2:nExpt,detFIMTrueInv_Unif,'-.r',...
    2:nExpt,detFIMTrueInv_Rand,'-.m',...
    'linewidth',2);

set(gca,"FontSize",16,'yscale','log')
legend({'|COV|','|FIM^-^1|', '|FIM^-^1| True','FIMopt','Uniform','Random'});
xlabel('Experiment Round')

%% Load DUSP1 Results
clear det*
load IterativeExperimentResults_DUSP1_simulated_FIMopt_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_FIM(i-1) = det(covLogMH{i}); % |COV| of MH samples
    detFIMInv_FIM(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]); % Predicted |FIM^-1| 
    detFIMTrueInv_FIM(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]); % True |FIM^-1|
    pars_FIM(i-1,:) = parametersFound{i}; % Parameter set found
end



load IterativeExperimentResults_DUSP1_simulated_intuition_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_int(i-1) = det(covLogMH{i}); % |COV| of MH samples
    detFIMInv_int(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]); % Predicted |FIM^-1| 
    detFIMTrueInv_int(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]); % True |FIM^-1|
    pars_int(i-1,:) = parametersFound{i}; % Parameter set found
end

load IterativeExperimentResults_DUSP1_simulated_uniform_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Unif(i-1) = det(covLogMH{i}); % |COV| of MH samples
    detFIMInv_Unif(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]); % Predicted |FIM^-1| 
    detFIMTrueInv_Unif(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]); % True |FIM^-1|
    pars_Unif(i-1,:) = parametersFound{i}; % Parameter set found
end

load IterativeExperimentResults_DUSP1_simulated_random_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Rand(i-1) = det(covLogMH{i}); % |COV| of MH samples
    detFIMInv_Rand(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]); % Predicted |FIM^-1| 
    detFIMTrueInv_Rand(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]); % True |FIM^-1|
    pars_Rand(i-1,:) = parametersFound{i}; % Parameter set found
end
%% Plot DUSP1 Results
figure(2); clf;
% Invisible plot lines for last legend entries
plot(2:nExpt,detCov_FIM,'k',...
    2:nExpt,detFIMInv_FIM,'--k',...
    2:nExpt,detFIMTrueInv_FIM,'-.k',...
    'linewidth',2);

hold on
% Plot |COV| of each Design
plot(2:nExpt,detCov_FIM,'b',...
    2:nExpt,detCov_int,'g',...
    2:nExpt,detCov_Unif,'r',...
    2:nExpt,detCov_Rand,'m',...
        'linewidth',2)

% Plot predicted FIM of each Design
plot(2:nExpt,detFIMInv_FIM,'--b',...
    2:nExpt,detFIMInv_int,'--g',...
    2:nExpt,detFIMInv_Unif,'--r',...
    2:nExpt,detFIMInv_Rand,'--m',...
        'linewidth',2)

% Plot true |FIM| of each Design 
plot(2:nExpt,detFIMTrueInv_FIM,'-.b', ...
    2:nExpt,detFIMTrueInv_int,'-.g',...
    2:nExpt,detFIMTrueInv_Unif,'-.r',...
    2:nExpt,detFIMTrueInv_Rand,'-.m',...
        'linewidth',2)


set(gca,"FontSize",16,'yscale','log')
legend('MH','Predicted','Exact')
% legend('FIM','Intuition','Uniform','Random')
legend({'|COV|','|FIM^-^1|', '|FIM^-^1| True','FIMopt','Intuition','Uniform','Random'});
xlabel('Experiment Round');

% Plot parameter set variation with each Experiment Round
figure(3); clf
plot(2:nExpt,pars_FIM,'b','linewidth',3); hold on
plot(2:nExpt,pars_Unif,'r','linewidth',3)
plot(2:nExpt,pars_Rand,'m','linewidth',3)
plot(2:nExpt,pars_int,'g','linewidth',3)
ylabel('Parameter Value')
xlabel('Experiment round')

simPars = [0.14,1,1,0.01]; % Simulated model Parameters
tmp = simPars;
% tmp = ([load('SGRS_model_v1.mat').SGRS_Model.parameters{1:4,2}]); %

for i=1:length(tmp)
    plot([2,nExpt],[tmp(i),tmp(i)],'k--','LineWidth',2)
end
legend({'FIMopt','','','',...
    'Uniform','','','',...
    'Random','','','','Intuition',...
    '','','',...
    'True Parameter Value'});
set(gca,"FontSize",16,'yscale','log')

% Plot Mean Fold Error with each Experiment Round
figure(4); clf
mean_fold_error_FIM = sum(abs(log(pars_FIM./tmp)),2);
mean_fold_error_Unif = sum(abs(log(pars_Unif./tmp)),2);
mean_fold_error_Rand = sum(abs(log(pars_Rand./tmp)),2);
mean_fold_error_int = sum(abs(log(pars_int./tmp)),2);

plot(2:nExpt,mean_fold_error_FIM,'linewidth',2); hold on
plot(2:nExpt,mean_fold_error_Unif,'linewidth',2);
plot(2:nExpt,mean_fold_error_Rand,'linewidth',2);
plot(2:nExpt,mean_fold_error_int,'linewidth',2);
xlabel('Experiment Round')
ylabel('Mean Fold Error')
legend('FIM','Uniform','Random','Intuition')
set(gca,"FontSize",16)


function predictedDetCov = predictCov(fimSet,parSet)
arguments
    fimSet
    parSet = [1:length(fimSet{1})];
end
detCovPreds = zeros(1,length(fimSet));
for i = 1:length(fimSet)
    detCovPreds(i) = det(fimSet{i}(parSet,parSet)^(-1));
end
predictedDetCov = median(detCovPreds);
end

