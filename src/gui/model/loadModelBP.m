function [app, event] = loadModelBP(app, event, array_file)
% function to load previously saved models into the GUI
if nargin<3
    % Creates a pop-up window to select the desired model
    [array_file,array_path] = uigetfile('*.mat','Select Model');
    addpath(array_path);
end

if isequal(array_file,0)
    return
else
    load(array_file, 'Mytable_model', 'myparameters');

    [~,struc_size] = size(myparameters); % Evaluates if there is more than one set of parameters saved to the model

    if struc_size > 1 % Having more than one set of parameters
        % Creates a pop-up figure with the ability to select the desired parameter set via a drop down menu
        myoptions = cell(struc_size,1);
        for i = 1:struc_size
            myoption = myparameters(i).name;
            myoptions(i) = cellstr(myoption);
        end
        app.ReactionsTabOutputs.paramVal = myparameters;
        f = figure('Position',[535 190 300 200]);

        mypopup = uicontrol(f,'Style','popup',...
            'String',myoptions,...
            'Position',[100 62 100 75],...
            'FontSize',12.5);
        set(mypopup,'Callback',{@loadParameterPopupCallback,app});
        uiwait(f)
        close(f);
        for i = 1:struc_size
            ind1(i) = strcmp({char(myparameters(i).name)},string(app.ReactionsTabOutputs.loadParams));
        end
        % ReactionsTabOutputs.inputs the selected parameters into the reactions tab
        index = find(ind1 == 1);

    else % Without more than one set of parameters, the values are just put into the reactions tab
        index=1;
    end
    
    % input citation in saved model into the citation box
    app.ReactionsTabOutputs.citations = {char(myparameters(index).citation)};
        
    app.ModelReactionTable.Data = Mytable_model;
    app.ReactionsTabOutputs.presetInputs = myparameters(index).inputs;
    app.ReactionsTabOutputs.presetParameters = myparameters(index).value;
    app.ReactionsTabOutputs.parameters = {};
    app.ReactionsTabOutputs.parameters(:,2) = app.ReactionsTabOutputs.presetParameters;
    try
        app.ReactionsTabOutputs.parameters(:,1) = myparameters(index).par_names;
    catch
        app.ReactionsTabOutputs.parameters(:,1)={''};
    end
    
    app.ModelParameterTable.Data = [myparameters(index).par_names,myparameters(index).value];
    
    % input initial consitions used for the ssa and fsp tabs
    app.ReactionsTabOutputs.initialCondition = char(myparameters(index).initialcond);

    % fills in the about the model box
    app.ModelAbout.Value=...
        {('About this model:');('');...
    ['Name of parameter set: ',char(myparameters(index).name)];('')};
    try
        app.ModelAbout.Value{end+1} = ['Model Description: ',char(myparameters(index).description)];
        app.ModelAbout.Value{end+1} = '';
    catch
    end
    try
        app.ModelAbout.Value{end+1} = ['Species Names: ',char(myparameters(index).species)];
        app.ModelAbout.Value{end+1} = '';
    catch
    end
    try
        app.ModelAbout.Value{end+1} = ['Curator: ',char(myparameters(index).curator)];
        app.ModelAbout.Value{end+1} = '';
    catch
    end
    try
        app.ModelAbout.Value{end+1} = ['Email: ',char(myparameters(index).email)];
        app.ModelAbout.Value{end+1} = '';
    catch
    end

end
end
