function [app, event] = loadModelBP(app, event, fileName)
arguments
    app = []
    event = []
    fileName = []
end
% function to load previously saved models into the GUI
if isempty(fileName)
    % Creates a pop-up window to select the desired model
    [array_file,array_path] = uigetfile('*.mat','Select Model');
    fileName = fullfile(array_path, array_file);
end

if isequal(fileName,0)
    disp('No file selected');
    return
else
    app.ModelFile.fileName = fileName;
    info = whos('-file', fileName);
    
    % Check if any variable is of class 'SSIT'
    if any(strcmp({info.class}, 'SSIT'))
        app.ChooseSSITModel.Visible = 'on';
        app.ChooseSSITModel_2.Visible = 'on';
        app.ChooseSSITModelLabel.Visible = 'on';
        app.ChooseSSITModelLabel_2.Visible = 'on';
        app.ChooseSSITModel.Items = {'Select Model',info(strcmp({info.class}, 'SSIT')).name};
        if length(app.ChooseSSITModel.Items)==2
            % Model = load(fileName,app.ChooseSSITModel.Items{1});
            % app.SSITModel = Model.(app.ChooseSSITModel.Items{1});
            app.ChooseSSITModel.Value = app.ChooseSSITModel.Items{2};
            app.ModelFile.modelName = app.ChooseSSITModel.Items{2};
            app = ChooseSSITModelValue(app);
            % updateAppFromSSIT(app);

            info = dir(app.ModelFile.fileName);
            app.FileModelLabel.Text = {['File: ',app.ModelFile.fileName];...
                ['Model: ',app.ModelFile.modelName];...
                ['Last Saved: ',info.date]};        
        end
        return
    end
    
    load(fileName, 'Mytable_model', 'myparameters');
    [~,struc_size] = size(myparameters); % Evaluates if there is more than one set of parameters saved to the model

    if struc_size > 1 % Having more than one set of parameters
        % % Creates a pop-up figure with the ability to select the desired parameter set via a drop down menu
        % myoptions = cell(struc_size,1);
        % for i = 1:struc_size
        %     myoption = myparameters(i).name;
        %     myoptions(i) = cellstr(myoption);
        % end
        % app.ReactionsTabOutputs.paramVal = myparameters;
        % f = figure('Position',[535 190 300 200]);
        % 
        % mypopup = uicontrol(f,'Style','popup',...
        %     'String',myoptions,...
        %     'Position',[100 62 100 75],...
        %     'FontSize',12.5);
        % set(mypopup,'Callback',{@loadParameterPopupCallback,app});
        % uiwait(f)
        % close(f);
        % for i = 1:struc_size
        %     ind1(i) = strcmp({char(myparameters(i).name)},string(app.ReactionsTabOutputs.loadParams));
        % end
        % % ReactionsTabOutputs.inputs the selected parameters into the reactions tab
        % index = find(ind1 == 1);
        % error('Support for multiple parameter sets has been removed');
        warning('Support for multiple parameter sets has been removed. Only using the first one.');
        index=1;
        % TODO - remove this section when all multi-par examples have been
        % converted.


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

    %% Generate and save SSIT model.
    app.SSITModel = SSIT('Empty');
    % Use default species names of x1, x2, ...
    
    % Detect species names from table (this version only supports x1, x2, ...)
    nSp = 0;
    allPlayers = append(Mytable_model{:,2:3});
    while contains(allPlayers,['x',num2str(nSp+1)])
        nSp = nSp+1;
        app.SSITModel.species{nSp,1} = ['x',num2str(nSp)];
    end

    % Build reaction network.
    nRxn = size(Mytable_model,1);
    app.SSITModel.propensityFunctions = cell(nRxn,1);
    app.SSITModel.stoichiometry = zeros(nSp,nRxn);
    for iRxn = 1:nRxn
        app.SSITModel.propensityFunctions{iRxn} = Mytable_model{iRxn,4};
        for iSp = 1:nSp
            if contains(Mytable_model{iRxn,2},app.SSITModel.species{nSp,1})
                J = strfind(Mytable_model{iRxn,2},app.SSITModel.species{nSp,1});
                K1 = strfind(Mytable_model{iRxn,2}(J:end),'(');
                K2 = strfind(Mytable_model{iRxn,2}(J:end),')');
                app.SSITModel.stoichiometry(iSp,iRxn) = -eval(Mytable_model{iRxn,2}(J+K1:J+K2-2));
            end
            if contains(Mytable_model{iRxn,3},app.SSITModel.species{nSp,1})
                J = strfind(Mytable_model{iRxn,3},app.SSITModel.species{nSp,1});
                K1 = strfind(Mytable_model{iRxn,3}(J:end),'(');
                K2 = strfind(Mytable_model{iRxn,3}(J:end),')');
                app.SSITModel.stoichiometry(iSp,iRxn) = app.SSITModel.stoichiometry(iSp,iRxn)+eval(Mytable_model{iRxn,3}(J+K1:J+K2-2));
            end
        end
    end
    
    app.SSITModel.parameters = [myparameters.par_names,myparameters.value];
    app.SSITModel.inputExpressions = [myparameters.input_names,myparameters.inputs];
    app.SSITModel.initialCondition = zeros(length(app.SSITModel.species),1);
    app.SSITModel.description = app.ModelAbout.Value;
    
    k = strfind(fileName,'.'); k=k(end);
    k1 = strfind(fileName,'/'); k1=k1(end);
    ModelName = fileName(k1+1:k-1);
    
    if strcmp(fileName(end-1:end),'.m')
        fileName = append(fileName(1:end-2),'.SSIT.mat');
    elseif strcmp(fileName(end-3:end),'.mat')
        fileName = append(fileName(1:end-4),'.SSIT.mat');
    end
    eval(append(ModelName,' = app.SSITModel;'));
    save(fileName,ModelName);
    updateAppFromSSIT(app);
    % updateModel(app);
end

end
