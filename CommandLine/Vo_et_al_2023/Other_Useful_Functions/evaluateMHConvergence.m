close all
% load TMPmh_5_667.mat
% load TMPmh_3_974.mat
% load TMPmh_5_BrianProcessingTS2_PDO1_550.mat
fName = 'TMPmh_1_Name_for_fitting_case';
load(fName)

% resave=0;
% for iPDO = 1:5
%     smpl = chainResults{iPDO}.mhSamples;
%     value = chainResults{iPDO}.mhValue;
%     smplDone = squeeze(smpl(value~=0,:));
%     valDone = squeeze(value(value~=0));
%     [~,j] = max(valDone)
%     if j~=1
%         ModelZero{iPDO}.parameters(ModelZero{iPDO}.fittingOptions.modelVarsToFit,2)=...
%             num2cell(exp(smplDone(j,:))');
%         resave=1;
%     end
% end
% if resave
%     save(fName,'FIMZero','FIMZeroLog','ModelZero','chainResults','covLog','fimResults','lTrue2FISH','lTrue2MCP','fimResultsOne')
% end
     
% smpl = chainResults{4}.mhSamples;
% value = chainResults{4}.mhValue;
smplDone = squeeze(smpl(value~=0,:));
valDone = squeeze(value(value~=0));

[~,j] = max(valDone)


% AAA = [exp(smplDone(j,:));...
% [ModelZero{1}.parameters{1:4,2}];...
% exp(smplDone(1,:))]

figure(1)
plot(valDone)
figure(2)
for i =1:size(smplDone,2)-1
    for j = i+1:size(smplDone,2)
        figure(2)
        subplot(size(smplDone,2)-1,size(smplDone,2)-1,(i-1)*(size(smplDone,2)-1)+j-1)
        plot(smplDone(:,j)/log(10),smplDone(:,i)/log(10),'ro-')
        figure(3)
        subplot(size(smplDone,2)-1,size(smplDone,2)-1,(i-1)*(size(smplDone,2)-1)+j-1)
        plot(exp(smplDone(:,j)),exp(smplDone(:,i)),'ro-')
    end
end

%% Effective Sample Size
figure(4); clf
for ipar = 1:size(smplDone,2)+1
    if ipar<=size(smplDone,2)
        ac = xcorr(smplDone(:,ipar)-mean(smplDone(:,ipar)),'normalized');
    else
        ac = xcorr(valDone-mean(valDone),'normalized');
    end
    ac = ac(size(smplDone,1):end);
    plot(ac,'LineWidth',3); hold on
    N = size(smplDone,1);
    tau = 1+2*sum((ac(2:N/100)));
    Neff(ipar) = N/tau;
end
Neff
legend('k12','k21','g','logL')

%%
