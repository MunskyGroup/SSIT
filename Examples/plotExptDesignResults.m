
% clear all
% close all
% parfor i=4:6
%     switch i
%         case 1
%             iterativeExperimentRunner('Poisson','simulated','FIMopt',10,1,i)
%         case 2
%             iterativeExperimentRunner('Poisson','simulated','FIMopt',10,1,i)
%         case 3
%             iterativeExperimentRunner('Poisson','simulated','uniform',10,1,i)
%         case 4            
%             iterativeExperimentRunner('DUSP1','simulated','FIMopt',10,1,i)
%         case 5
%             iterativeExperimentRunner('DUSP1','simulated','random',10,1,i)
%         case 6
%             iterativeExperimentRunner('DUSP1','simulated','uniform',10,1,i)
%     end
% end


%% Poisson Results
% clear det*
% load IterativeExperimentResults_Poisson_simulated_FIMopt_2
% nExpt = length(covLogMH);
% for i = 2:nExpt
%     detCov_FIM(i-1) = det(covLogMH{i});
%     detFIMInv_FIM(i-1) = predictCov(FIMcurrentExptSaved{i});
%     detFIMTrueInv_FIM(i-1) = predictCov(FIMcurrentExptTrueSaved{i});
% end
% 
% load IterativeExperimentResults_Poisson_simulated_uniform_2
% nExpt = length(covLogMH);
% for i = 2:nExpt
%     detCov_Unif(i-1) = det(covLogMH{i});
%     detFIMInv_Unif(i-1) = predictCov(FIMcurrentExptSaved{i});
%     detFIMTrueInv_Unif(i-1) = predictCov(FIMcurrentExptTrueSaved{i});
% end
% 
% load IterativeExperimentResults_Poisson_simulated_random_2
% nExpt = length(covLogMH);
% for i = 2:nExpt
%     detCov_Rand(i-1) = det(covLogMH{i});
%     detFIMInv_Rand(i-1) = predictCov(FIMcurrentExptSaved{i});
%     detFIMTrueInv_Rand(i-1) = predictCov(FIMcurrentExptTrueSaved{i});
% end
% 
% figure(1); clf;
% plot(2:nExpt,detCov_FIM,'b',...
%     2:nExpt,detCov_Unif,'r',...
%     2:nExpt,detCov_Rand,'m',...
%     'linewidth',2)
% 
% hold on
% plot(2:nExpt,detFIMInv_FIM,'--b',...
%     2:nExpt,detFIMInv_Unif,'--r',...
%     2:nExpt,detFIMInv_Rand,'--m',...
%     'linewidth',2)
% 
% plot(2:nExpt,detFIMTrueInv_FIM,'-.b',...
%     2:nExpt,detFIMTrueInv_Unif,'-.r',...
%     2:nExpt,detFIMTrueInv_Rand,'-.m',...
%     'linewidth',2)
% 
% set(gca,"FontSize",16,'yscale','log')
% legend('FIM','Uniform','Random')

%% DUSP1 Results
clear det*
% load IterativeExperimentResults_DUSP1_simulated_FIMopt_4
load IterativeExperimentResults_DUSP1_real_FIMopt_4

nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_FIM(i-1) = det(covLogMH{i});
    detFIMInv_FIM(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]);
    detFIMTrueInv_FIM(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]);
    pars_FIM(i-1,:) = parametersFound{i};
end

% load IterativeExperimentResults_DUSP1_simulated_intuition_4
load IterativeExperimentResults_DUSP1_real_intuition_4

nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_int(i-1) = det(covLogMH{i});
    detFIMInv_int(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]);
    detFIMTrueInv_int(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]);
    pars_int(i-1,:) = parametersFound{i};
end

% load IterativeExperimentResults_DUSP1_simulated_uniform_4
load IterativeExperimentResults_DUSP1_real_uniform_4
nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Unif(i-1) = det(covLogMH{i});
    detFIMInv_Unif(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]);
    detFIMTrueInv_Unif(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]);
    pars_Unif(i-1,:) = parametersFound{i};
end

% load IterativeExperimentResults_DUSP1_simulated_random_4
load IterativeExperimentResults_DUSP1_real_random_4

nExpt = length(covLogMH);
for i = 2:nExpt
    detCov_Rand(i-1) = det(covLogMH{i});
    detFIMInv_Rand(i-1) = predictCov(FIMcurrentExptSaved{i},[1:4]);
    detFIMTrueInv_Rand(i-1) = predictCov(FIMcurrentExptTrueSaved{i},[1:4]);
    pars_Rand(i-1,:) = parametersFound{i};
end

figure(2); clf;
plot(2:nExpt,detCov_FIM,'b',...
    2:nExpt,detCov_Unif,'r',...
    2:nExpt,detCov_Rand,'m',...
    2:nExpt,detCov_int,'g',...
    'linewidth',2)

hold on
plot(2:nExpt,detFIMInv_FIM,'--b',...
    2:nExpt,detFIMInv_Unif,'--r',...
    2:nExpt,detFIMInv_Rand,'--m',...
    2:nExpt,detFIMInv_int,'--g',...
    'linewidth',2)

plot(2:nExpt,detFIMTrueInv_FIM,'-.b',...
    2:nExpt,detFIMTrueInv_Unif,'-.r',...
    2:nExpt,detFIMTrueInv_Rand,'-.m',...
    2:nExpt,detFIMTrueInv_int,'-.g',...
    'linewidth',2)

set(gca,"FontSize",16,'yscale','log')
legend('MH','Predicted','Exact')
legend('FIM','Uniform','Random','Intuition')
xlabel('Experiment round');

figure(3); clf
plot(2:nExpt,pars_FIM,'b','linewidth',3); hold on
plot(2:nExpt,pars_Unif,'r','linewidth',3)
plot(2:nExpt,pars_Rand,'m','linewidth',3)
plot(2:nExpt,pars_int,'g','linewidth',3)

figure(4); clf
mean_fold_error_FIM = sum(abs(log(pars_FIM./tmp)),2);
mean_fold_error_Unif = sum(abs(log(pars_Unif./tmp)),2);
mean_fold_error_Rand = sum(abs(log(pars_Rand./tmp)),2);
mean_fold_error_int = sum(abs(log(pars_int./tmp)),2);

plot(2:nExpt,mean_fold_error_FIM); hold on
plot(2:nExpt,mean_fold_error_Unif);
plot(2:nExpt,mean_fold_error_Rand);
plot(2:nExpt,mean_fold_error_int);


tmp = ([load('SGRS_model_v1.mat').SGRS_Model.parameters{1:4,2}]); %
for i=1:length(tmp)
    plot([2,nExpt],[tmp(i),tmp(i)],'k--','LineWidth',2)
end
set(gca,"FontSize",16,'yscale','log')


% legend(legend(['MH','Predicted','Exact';'FIM','Uniform','Random']);
% %%
% return
% %%
% load IterativeExperimentResults_DUSP1_real_FIMopt_4
% fim_ = predictCov(FIMcurrentExptTrueSaved{2},[1:4])
% 
% load IterativeExperimentResults_DUSP1_real_uniform_4.mat
% uni_ = predictCov(FIMcurrentExptTrueSaved{2},[1:4])
% 
% load IterativeExperimentResults_DUSP1_real_intuition_4.mat
% int_ = predictCov(FIMcurrentExptTrueSaved{2},[1:4])
% 
% load IterativeExperimentResults_DUSP1_real_random.mat
% ran_ = predictCov(FIMcurrentExptTrueSaved{2},[1:4])


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

% plot(2:nExpt,detFIMTrueInv_FIM,'k',2:nExpt,detFIMTrueInv_FIM,'--k',2:nExpt,detFIMTrueInv_FIM,'-.k','linewidth',2)
% legend('|COV|','|FIM^-1|','|FIM^-1| True')
