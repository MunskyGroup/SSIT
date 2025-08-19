function app = clearDataConstraints(app)
for iConstr = 1:6
    if iConstr>1
        app.(['AndOr',num2str(iConstr)]).Value = '-';
    end
    app.(['DataConstrText',num2str(iConstr)]).Value = '-';

    if ~isempty(app.DataLoadingAndFittingTabOutputs.columns)
        app.(['DataConstrChoice',num2str(iConstr)]).Items = ['-';app.DataLoadingAndFittingTabOutputs.columns'];
    else
        app.(['DataConstrChoice',num2str(iConstr)]).Items = {'-'};
    end
    app.(['DataConstrChoice',num2str(iConstr)]).Value = '-';
end
app = CreateDataConstaintOptions(app);