function SetPriors(app)
app.DataLoadingAndFittingTabOutputs.priorOptions = ssit.parest.propsStorage;
clear TMP F
switch app.PriorTypeDropDown.Value
    case 'None'
        app.SpecifyPrioronParametersButton.Enable = 0;
        TMP=[];
    case 'Normal'
        app.SpecifyPrioronParametersButton.Enable = 1;
        J = find(strcmp(app.fit_parameters_table.Data(:,3),'y'));
        for i=1:length(J)
            TMP(i*2-1,:) = {['MEAN_',app.fit_parameters_table.Data{i,1}],'100'};
            TMP(i*2,:) = {['STDV_',app.fit_parameters_table.Data{i,1}],'10'};
        end
        if ~isempty(app.pdo_parameters_table.Data)
            J = find(strcmp(app.pdo_parameters_table.Data(:,3),'y'));
            for k=1:length(J)
                TMP((i+k)*2-1,:) = {['MEAN_',app.pdo_parameters_table.Data{k,1}],'0.5'};
                TMP((i+k)*2,:) = {['STDV_',app.pdo_parameters_table.Data{k,1}],'0.1'};
            end
        end
    case 'LogNormal'
        app.SpecifyPrioronParametersButton.Enable = 1;
        J = find(strcmp(app.fit_parameters_table.Data(:,3),'y'));
        for i=1:length(J)
            TMP(i*2-1,:) = {['LogMEAN_',app.fit_parameters_table.Data{i,1}],'0'};
            TMP(i*2,:) = {['LogSTDV_',app.fit_parameters_table.Data{i,1}],'1'};
        end
        if ~isempty(app.pdo_parameters_table.Data)
            J = find(strcmp(app.pdo_parameters_table.Data(:,3),'y'));
            for k=1:length(J)
                TMP((i+k)*2-1,:) = {['LogMEAN_',app.pdo_parameters_table.Data{k,1}],'0'};
                TMP((i+k)*2,:) = {['LogSTDV_',app.pdo_parameters_table.Data{k,1}],'1'};
            end
        end
end

F = struct;
for i = 1:size(TMP,1)
    F.(TMP{i,1}) = TMP{i,2};
end

app.DataLoadingAndFittingTabOutputs.priorOptions.props = F;
end