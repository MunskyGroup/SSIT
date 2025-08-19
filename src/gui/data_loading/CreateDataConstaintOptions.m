function app = CreateDataConstaintOptions(app)
% This function updates options for data constaints based on type of data
% chosen for the constraint.
for iConstr = 1:6
    J = find(strcmp(app.DataLoadingAndFittingTabOutputs.columns,app.(['DataConstrChoice',num2str(iConstr)]).Value));
    if isempty(J)
        app.(['DataLogical',num2str(iConstr)]).Items = {'NA'};
    else
        switch app.DataLoadingAndFittingTabOutputs.columnTypes{J}
            case 'numeric'
                app.(['DataLogical',num2str(iConstr)]).Items = {'==';'~=';'>';'>=';'<';'<='};
            case 'string'
                app.(['DataLogical',num2str(iConstr)]).Items = {'is =';'is not ='};
        end
    end
end