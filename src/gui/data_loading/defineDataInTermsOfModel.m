function defineDataInTermsOfModel(app)
%% This function selects the data from the loaded file and assigns the the 
% data in the file as the appropiate name depending on the selected boxes.
% The user selects the checkboxes and drop down value to assign the data to
% x1, x2, or x3.
%% Apply filter to select data from user specified conditions
if app.FilteronConditionCheckBox.Value
    J_col = find(strcmp(app.cond_DropDown.Items,app.cond_DropDown.Value))';
    cond = app.DataLoadingAndFittingTabOutputs.dataTable(:,J_col);
    if iscell(cond)  % are data conditions in cells
        if iscellstr(cond(1))  % Do these cells contain strings
            cond_to_keep = 'str';
        elseif isnumeric(cond{1})  % These cells contain numbers
            cond_to_keep = 'numcell';
        end
        % Conditions are not string
    elseif isnumeric(cond(1))
        cond_to_keep = 'num';
    end
    
    switch cond_to_keep
        case 'str'
            J_row = strcmp(app.con_DropDown_2.Value,cond);
        case 'num'
            J_row = str2num(app.con_DropDown_2.Value)==cond;
        case 'numcell'
            J_row = str2num(app.con_DropDown_2.Value)==cell2mat(cond);
    end
    app.DataLoadingAndFittingTabOutputs.filteredDataMatrix = app.DataLoadingAndFittingTabOutputs.dataTable(J_row,:);
else
    app.DataLoadingAndFittingTabOutputs.filteredDataMatrix = app.DataLoadingAndFittingTabOutputs.dataTable;
end

%% Determine which columns of data correspond to which species.
app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(1) = find(contains(app.ParEstX1DropDown.Items,app.ParEstX1DropDown.Value));
if app.ParEstX1CheckBox.Value==0||strcmp(app.ParEstX1DropDown.Value,'ignore')
    app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(1)=0;
end
app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(2) = find(contains(app.ParEstX1DropDown.Items,app.ParEstX2DropDown.Value));
if app.ParEstX2CheckBox.Value==0||strcmp(app.ParEstX2DropDown.Value,'ignore')
    app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(2)=0;
end
app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(3) = find(contains(app.ParEstX1DropDown.Items,app.ParEstX3DropDown.Value));
if app.ParEstX3CheckBox.Value==0||strcmp(app.ParEstX3DropDown.Value,'ignore')
    app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(3)=0;
end

%% Separate Data according to different measurement times.
list_of_data_times = cell2mat(app.DataLoadingAndFittingTabOutputs.filteredDataMatrix(:,app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_time_index));
% whos list_of_data_times
% size(app.DataLoadingAndFittingTabOutputs.filteredDataMatrix)
% size(app.DataLoadingAndFittingTabOutputs.dataTable)
app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_data_times = {};
for i = 1:length(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)
    app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_data_times{i} = find(list_of_data_times == app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times(i));
end
         
%% Separate Data into Tensor Formulation.
%%      Initialize tensor matrix to appropriate size
tensor_size = ones(1,4);
tensor_size(1) = length(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times);
for i = 1:3
    if app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(i)>0
        try
            tensor_size(i+1) = max(cell2mat(app.DataLoadingAndFittingTabOutputs.filteredDataMatrix(:,app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(i))))+1;
        catch
            app.DataLoadingAndFittingTabOutputs.filteredDataMatrix(:,app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(i))
            tensor_size(i+1) = max((app.DataLoadingAndFittingTabOutputs.filteredDataMatrix(:,app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(i))))+1;
        end
    else
        tensor_size(i+1)=1;
    end
end
app.DataLoadingAndFittingTabOutputs.dataTensor = zeros(tensor_size);
% %%      Sort data counts into tensor
% X = app.DataLoadingAndFittingTabOutputs.filteredDataMatrix;
% J = app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species;
% X = X(:,[true,J~=0]);
% for it = 1:length(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)
%     XX = X(app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_data_times{it},2:end);
%     Y = zeros(size(XX));
%     Y(:) = [XX{:}];
%     [C,~,IC] = unique(Y,'rows','stable');
%     for ic = 1:size(C,1)
%         app.DataLoadingAndFittingTabOutputs.dataTensor(it,C(ic)+1) = sum(IC==ic);
%     end
% end
% 
%%
% for it = 1:length(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)
%     for icc = 1:length(app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_data_times{it})
%         ic = app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_data_times{it}(icc);
%         J = ones(1,3);
%         for is = 1:3
%             if app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(is)>0
%                 J(is) = cell2mat(app.DataLoadingAndFittingTabOutputs.filteredDataMatrix(ic,app.DataLoadingAndFittingTabOutputs.fittingOptions.indices_of_species(is)))+1;
%             end
%         end
%         app.DataLoadingAndFittingTabOutputs.dataTensor(it,J(1),J(2),J(3)) = app.DataLoadingAndFittingTabOutputs.dataTensor(it,J(1),J(2),J(3))+1;
%     end
% end
%%      Sort data counts into tensor
% This works, but is very slow.
tmpv = app.DataLoadingAndFittingTabOutputs.fittingOptions;
tmplog = find(tmpv.indices_of_species>0);
for it = 1:length(tmpv.fit_times)
    for icc = 1:length(tmpv.indices_of_data_times{it})
        ic =tmpv.indices_of_data_times{it}(icc);
        J = ones(1,3);
%         for is = tmplog%1:3
%             if tmpv.indices_of_species(is)>0
                J(tmplog) = cell2mat(app.DataLoadingAndFittingTabOutputs.filteredDataMatrix(ic,tmpv.indices_of_species(tmplog)))+1;
%             end
%         end
        app.DataLoadingAndFittingTabOutputs.dataTensor(it,J(1),J(2),J(3)) = app.DataLoadingAndFittingTabOutputs.dataTensor(it,J(1),J(2),J(3))+1;
    end
end
%%
size(app.DataLoadingAndFittingTabOutputs.dataTensor)