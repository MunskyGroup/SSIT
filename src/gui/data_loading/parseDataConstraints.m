function parseDataConstraints(app)

nSp = length(app.SSITModel.species);
app.DataLoadingAndFittingTabOutputs.linking = {};
for iSp = 1:nSp
    if ~strcmp(app.(['DataSpecies',num2str(iSp)]).Value,'none')
        app.DataLoadingAndFittingTabOutputs.linking(end+1,1:2) = {app.SSITModel.species{iSp},app.(['DataSpecies',num2str(iSp)]).Value};
    end
end

constraints='';
if ~strcmp(app.DataConstrChoice1.Value,'-')&&~strcmp(app.DataConstrText1.Value,'-')
    addConstr = true;
    switch app.DataLogical1.Value
        case 'is ='
            constraints = [constraints,'(strcmp(TAB.',app.DataConstrChoice1.Value,',',app.DataConstrText1.Value,')'];
        case 'is not ='
            constraints = [constraints,'(~strcmp(TAB.',app.DataConstrChoice1.Value,',',app.DataConstrText1.Value,')'];
        case {'==','~=','>','<','>=','<='}
            constraints = [constraints,'(TAB.',app.DataConstrChoice1.Value,app.DataLogical1.Value,app.DataConstrText1.Value,')'];
    end
    iConstr = 2;
    while addConstr
        if ~strcmp(app.(['AndOr',num2str(iConstr)]).Value,'-')&&~strcmp(app.(['DataConstrChoice',num2str(iConstr)]).Value,'-')&&~strcmp(app.(['DataConstrText',num2str(iConstr)]).Value,'-')
            switch app.(['AndOr',num2str(iConstr)]).Value
                case 'and'
                    conx = '&';
                case 'or'
                    conx = '|';              
            end
            constraints = [constraints,conx];
            switch app.(['DataLogical',num2str(iConstr)]).Value
                case 'is ='
                    constraints = [constraints,'(strcmp(TAB.',app.(['DataConstrChoice',num2str(iConstr)]).Value,',',app.(['DataConstrText',num2str(iConstr)]).Value,')'];
                case 'is not ='
                    constraints = [constraints,'(~strcmp(TAB.',app.(['DataConstrChoice',num2str(iConstr)]).Value,',',app.(['DataConstrText',num2str(iConstr)]).Value,')'];
                case {'>','<','>=','<='}
                    constraints = [constraints,'(TAB.',app.(['DataConstrChoice',num2str(iConstr)]).Value,app.(['DataLogical',num2str(iConstr)]).Value,app.(['DataConstrText',num2str(iConstr)]).Value,')'];
            end
            iConstr = iConstr+1;
        else
            addConstr=false;
        end
    end


end
if isempty(constraints)
    app.DataLoadingAndFittingTabOutputs.constraints = {};
    app.SSITModel = app.SSITModel.loadData(...
        app.DataLoadingAndFittingTabOutputs.dataFileName,...
        app.DataLoadingAndFittingTabOutputs.linking);
else
    app.DataLoadingAndFittingTabOutputs.constraints = {[],[],constraints};
    app.SSITModel = app.SSITModel.loadData(...
        app.DataLoadingAndFittingTabOutputs.dataFileName,...
        app.DataLoadingAndFittingTabOutputs.linking,...
        app.DataLoadingAndFittingTabOutputs.constraints);
end

app.NumberafterConstraintsLabel.Text = ['Number after Constraints: ',num2str(sum(app.SSITModel.dataSet.nCells))];
end