function [] = selectPresetExampleDropDown(app)
% This function autopopulates the preset models for the drop down menu
% within the reaction tab.
value = app.ModelUsePresetExampleTypeDropDown.Value;   % Evaluates the Selection of the Model type
% Finds the folder containing the preset examples within the correct OS
MatlabFiles = what(append('Models',filesep,value));

v2 = [MatlabFiles(end).m; 
    MatlabFiles(end).mat];

if length(MatlabFiles)>1
    whereAmI = pwd;
    for i = 1:length(MatlabFiles)
        J = min(length(whereAmI),length(MatlabFiles(end).path));
        K = max(length(whereAmI),length(MatlabFiles(end).path));
        matchNum(i) = sum(whereAmI(1:J)==MatlabFiles(end).path(1:J))/K;
    end
    [~,J] = max(matchNum);
    addpath(MatlabFiles(J).path)
    disp(['Warning-Multple SSIT directorys detected in path.  Using path = ',MatlabFiles(J).path])
end

app.ModelDropDown.Items = v2; % Fills the second dropdown menu to select the models within the category

end