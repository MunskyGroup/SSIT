function [] = selectPresetExampleDropDown(app)
% This function autopopulates the preset models for the drop down menu
% within the reaction tab.
value = app.ModelUsePresetExampleTypeDropDown.Value;   % Evaluates the Selection of the Model type
% Finds the folder containing the preset examples within the correct OS
% if ispc
%     MatlabFiles = what(['Models\',value]);    
% elseif isunix
    MatlabFiles = what(['Models/',value]);    
% else    
%     disp('Platform not supported');
% end

v2 = [MatlabFiles(end).m; 
    MatlabFiles(end).mat];

if length(MatlabFiles)>1
    disp(['Warning-Multple SSIT directorys detected in path.  Using path = ',MatlabFiles(end).path])
end

app.ModelDropDown.Items = v2; % Fills the second dropdown menu to select the models within the category

end