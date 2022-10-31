function randomGreedyFit(app)
for kit = 1:1000
    %% Uses older parameter values and FSP outputs to create the parameter 
    % values in the table.
    disp(['Running parameter try ',num2str(kit)]);
    Xold = cell2mat(app.tab_pars_Dan.Data(:,2));
    Pold = app.FspTabOutputs.solutions;
    
    X = Xold.*(1 + 0.1*randn(size(Xold)).*(rand(size(Xold))>0.5));
    app.tab_pars_Dan.Data(:,2) = num2cell(X);
    
    app.ModelParameterTable.Data=app.tab_pars_Dan.Data;
    
    if isempty(app.FspConstraintTable.Data)
        DefaultButtonPushed(app)
    end
    runFsp(app)
    updatePlotsinFsp(app)
    %%
    if strcmp(app.SpeciesDropDown.Value,'ILN1')
        K_spec = 1;
    else
        K_spec = 2;
    end
    
    T_array = eval(app.FspPrintTimesField.Value);
    times = 10+[0,0.5,1,2,4];
    tms = {'0hr','30mins','1hr','2hr','4hr'};
    for i=1:5
        str = ['outputs_',tms{i},'_',app.DrugDropDown.Value];
        load('Data/Dan_and_Jim/FinalDataAll.mat',str);
        D =  eval(str);
        
        [~,j] = min(abs(T_array-times(i)));
        
        Z = squeeze(app.FspTabOutputs.solutions(j,:,:,:));
        P = squeeze(sum(sum(Z,1),2));
        
        P = max(P,1e-6);
        D = D(K_spec,:)+1;
        D(D>length(P))=length(P);
        if length(P)<max(D)+1
            P(length(P)+1:max(D)+1) = 1e-6;
        end
        LogLk(i) = sum(log(P(D)));
    end
    if sum(LogLk)>app.DataLoadingAndFittingTabOutputs.J_LogLk
        Update_Dans_Plots(app);
        saveModelBP(app,event,app.NameForModelSavingEditField.Value)
    else
        
        app.tab_pars_Dan.Data(:,2) = num2cell(Xold);
        app.FspTabOutputs.solutions = Pold;
        
        app.ModelParameterTable.Data=app.tab_pars_Dan.Data;
%         UITableCellEdit(app)
    end
end
end