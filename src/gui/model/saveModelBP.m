function [app, event, par_name] = saveModelBP(app,event,override_name)
% Function saves the edits created to a Preset model or your own model

% defines questions in popup window when saving model
title = 'Save your model';
        
dims2 = [1 200];

prompt = {'Delete Previous Models? (yes or no)','no';...
    'Name of new parameter set (No spaces)','Set_1';...
    'Model Reference Article','not entered';...
    'Model Description','not entered';...
    'URL','not entered';...
    'Species Names','(x1=??,x2=??,x3=??)';...
    'Initial Conditions','[0,0,0]';...
    'SSIT Model Curator','not entered';...
    'email','not entered'};

if nargin<3
    % If saved previously, a search box will open to find the model as it was previously saved
    [array_file,directory] = uiputfile('*.mat','Save Model as');
    modelpath = string([directory,array_file]);                 % Creates a path to the same filename and folder as previously selected

    try
        load(modelpath,'myparameters');            % Opens the array file with the variables Mytable_model and myparameters
        [~,max_ind] = size(myparameters);                           % Identifies the final index of myparameters
        index = max_ind+1;                                          % Adds another field to myparameters in order to add new parameters to it
        definput2 = prompt(:,2);
        definput2(2) = {['Set_',num2str(index)]};
        answer = inputdlg(prompt(:,1),title,dims2,definput2);    
    catch
        index=1;
        answer= inputdlg(prompt(2:end,1),title,dims2,prompt(2:end,2));  
        answer=['yes';answer];
    end
else
    if exist(override_name,'file')
        load(override_name,'Mytable_model','myparameters');         % Opens the array file with the variables Mytable_model and myparameters
        [~,max_ind] = size(myparameters);                           % Identifies the final index of myparameters
        index = max_ind+1;                                          % Adds another field to myparameters in order to add new parameters to it
    else
        index=1;
    end
    myparameters(index).likelihood = app.DataLoadingAndFittingTabOutputs.J_LogLk;
    modelpath=override_name;
end
Mytable_model = app.ModelReactionTable.Data;                      % Sets the model as the reactions table of data
myparameters(index).value = app.ModelParameterTable.Data(:,2);          % Adds the values of the parameters to the new index of the structure
myparameters(index).inputs = app.ModelInputTable.Data(:,2);  % Adds the inputs of the new set of parameters to the new index of the structure
myparameters(index).par_names = app.ModelParameterTable.Data(:,1);
myparameters(index).input_names = app.ModelInputTable.Data(:,1);

myparameters(index).name = string(answer(2));               % Adds the name of the parameters to the new index of the structure
myparameters(index).citation = answer(3); 
myparameters(index).description = answer(4); 
myparameters(index).url = answer(5); 
myparameters(index).species = answer(6);
myparameters(index).initialcond = answer(7);
myparameters(index).curator = answer(8); 
myparameters(index).email = answer(9); 
answer

save(modelpath,'Mytable_model','myparameters')              % Re-saves the model with the new parameters

par_name = myparameters(index).name;

end
