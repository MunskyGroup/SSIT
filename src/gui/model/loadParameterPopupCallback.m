function loadParameterPopupCallback (hobject, event, app)
%% This function loads the parameter set for the model. 
% When a model is selected and there is more then one parameter set saved,
% the user will be prompted to select the parameter set from a drop down
% menu.

%     keyboard
    items = get(hobject,'String');
    index_selected = get(hobject,'Value');
    item_selected = items{index_selected};

    app.ReactionsTabOutputs.loadParams = item_selected;
    
    uiresume
   
end
