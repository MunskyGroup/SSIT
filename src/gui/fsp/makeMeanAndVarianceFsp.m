function makeMeanAndVarianceFsp(app)
% This function creates graphs of the different marginals at the desired
% time point in the FSP 3 Species GUI. This will use the time selected from
% the time slider as the time point from each of the different combinations
% based on the species selected.

T_array = eval(app.FspPrintTimesField.Value);
for it = 1:length(T_array)
    Nd = app.FspTabOutputs.solutions{it}.p.dim;
    for i=1:Nd
        if Nd==1
            mdist{1} = double(app.FspTabOutputs.solutions{it}.p.data);
        else
            INDS = setdiff([1:Nd],i);
            mdist{i} = double(app.FspTabOutputs.solutions{it}.p.sumOver(INDS).data);
        end
    end
    for j=1:Nd
        mns(it,j) = [0:length(mdist{j})-1]*mdist{j};
        mns2(it,j) = [0:length(mdist{j})-1].^2*mdist{j};
    end
end
vars = mns2-mns.^2;

figure();
Nd = app.FspTabOutputs.solutions{1}.p.dim;
for iSpecies = 1:Nd
    Plts_to_make(iSpecies) = max(contains(app.SpeciestoShowListBoxMeans.Value,app.SSITModel.species{iSpecies}));
end
cols = {'b','r','g','m','c','y','k'};
cols2 = {'--m','--c','--k','b--','r--','g--','m--','c--','y--','k--'};
LG = {};
for iplt=1:Nd
    if Plts_to_make(iplt)
        if app.FspMeanVarShowVarianceCheckBox.Value == 1
            shadedErrorBar(T_array,mns(:,iplt),sqrt(vars(:,iplt)),'lineprops',cols{iplt});hold on;
            LG{end+1} = [app.SSITModel.species{iplt},' \pm std'];
        else
            plot(T_array,mns(:,iplt),cols{iplt});hold on;
            LG{end+1} = [char(app.NameTable.Data(iplt,2))];
        end
        if app.FspMeanVarShowOdeCheckBox.Value == 1
            if isempty(app.FspTabOutputs.odeSolutions)
                runOde(app);
            end
            plot(app.FspTabOutputs.tOde,app.FspTabOutputs.odeSolutions(:,iplt),cols2{iplt},'linewidth',2); hold on;
            plot(T_array,mns(:,iplt),cols{iplt});
            LG{end+1} = [append(app.SSITModel.species{iplt},' ode')];
            LG{end+1} = app.SSITModel.species{iplt};
        end
    end
end
set(gca,'fontsize',20)
title('Averages Response vs Time');
xlabel('Time');
ylabel('Average Response');
legend(LG)

end
