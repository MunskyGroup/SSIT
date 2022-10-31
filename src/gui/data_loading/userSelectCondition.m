function userSelectCondition(app)
%% This function allows the user to choose what the selected value will 
% be represented as in the drop down menu.
if app.FilteronConditionCheckBox.Value
    app.cond_DropDown.Items = app.ParEstX1DropDown.Items(1:end-1);
    J_col = strcmp(app.cond_DropDown.Items,app.cond_DropDown.Value);
    try
        choices = (unique(app.DataLoadingAndFittingTabOutputs.dataTable(:,J_col)));
    catch
        choices = unique(cell2mat(app.DataLoadingAndFittingTabOutputs.dataTable(:,J_col)));
    end
    
    try
        for i=1:length(choices)
            choice_cell{i} = num2str(choices(i));
        end
    catch
        choice_cell= choices;
    end
    
    app.con_DropDown_2.Items = choice_cell;
    defineDataInTermsOfModel(app)
else
    app.cond_DropDown.Items = {'keep all'};
    app.cond_DropDown.Value = {'keep all'};
    app.con_DropDown_2.Items = {'empty'};
    app.con_DropDown_2.Value = {'empty'};
    defineDataInTermsOfModel(app)
end
   
