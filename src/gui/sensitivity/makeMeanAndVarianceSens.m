function makeMeanAndVarianceSens(app)
% this fucntion takes the derivitive of the mean and variance and displays
% it for each species.

% Find the time index to plot
T_array = eval(app.SensPrintTimesEditField.Value);

%% Compute the marginal distributions and Find mean and variance of marginals
solutionFormat = app.SensFspTabOutputs.solutions.format;

for it = 1:length(T_array)
    if (strcmp(solutionFormat, 'forward'))
        mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{it}.states, ...
            app.SensFspTabOutputs.solutions.data{it}.p,false);
        %     sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states,...
        %         app.SensFspTabOutputs.solutions.data{j}.dp(:, ipar),false);
    else
        mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.ps{it}.states, ...
            app.SensFspTabOutputs.solutions.data.ps{it}.p,false);
        %     sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.dps{j, ipar}.states,...
        %         app.SensFspTabOutputs.solutions.data.dps{j, ipar}.p,false);
    end
    for j=1:3
        mns(it,j) = [0:length(mdist{j})-1]*mdist{j};
        mns2(it,j) = [0:length(mdist{j})-1].^2*mdist{j};
    end
end
Var = mns2 - mns.^2;

%% Take derivative of the mean and variance and plot
for k = 1:3
    mDiff(:,k) = gradient(mns(:,k));
    varDiff(:,k) = gradient(Var(:,k));
end
% keyboard
% make plots

figure();
Plts_to_make = [app.SensMeanAndVarX1CheckBox.Value,app.SensMeanAndVarX2CheckBox.Value,app.SensMeanAndVarX3CheckBox.Value];
cols = ['b','r','g'];
LG = {};
for iplt=1:3
    if Plts_to_make(iplt)
        if app.SensShowVarianceCheckBox.Value == 1
            shadedErrorBar(T_array,mDiff(:,iplt),sqrt(varDiff(:,iplt)),'lineprops',cols(iplt));hold on;
            LG{end+1} = [char(app.NameTable.Data(iplt,2)),' \pm std'];
        else
            plot(T_array,mDiff(:,iplt),cols(iplt));hold on;
            LG{end+1} = [char(app.NameTable.Data(iplt,2))];
        end
    end
end
set(gca,'fontsize',20)
title('Averages Response vs Time');
xlabel('Time');
ylabel('Average Response');
legend(LG)