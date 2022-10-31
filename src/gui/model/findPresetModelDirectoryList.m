function [event] = findPresetModelDirectoryList(app)
% This function occurs when the app starts up and populates the first drop 
% down menu with the elements from the directory Models

    v = dir('Models');                                          % Opens the folders found in the models folder to load the pre-set models
    % Remove non-directories for list.
    I = ones(1,length(v));
    for i=1:length(v)
        if v(i).name(1) == '.'
            I(i)=0;
        end
    end
    v = v(I==1);                                                % The first two indeces are empty, so the third to the end folder are the only titles loaded

    app.ModelUsePresetExampleTypeDropDown.Items = {v.name};     % This sets the items of the first drop down menu of the preset examples
    event = app.ModelUsePresetExampleTypeDropDown.Items(1);     % This selects and populates the first element of the drop down on starting
end        