%% Start app and create links to call back functions
close all force
clear all
cd ..
A = SSIT_GUI;  % Start up the app.


mod_change_d = get(eval(['A.ModelUsePresetExampleTypeDropDown']),'ValueChangedFcn');
mod_change = get(eval(['A.ModelDropDown']),'ValueChangedFcn');
run_SSA = get(eval('A.SsaRunButton'),'ButtonPushedFcn'); 
data_save = get(eval('A.SsaExportExcelButton'),'ButtonPushedFcn');
load_data = get(eval('A.ParEstLoadDataButton'),'ButtonPushedFcn');
run_fit = get(eval('A.FitModelButton'),'ValueChangedFcn');

%% 1) Make list of models to test
Models = {'epidemiology','sir_model_r0.mat';...
    'gene_expression','Gene_Toggle_Model.m';...
    'social','opinion_formation.mat';...
    'gene_expression','poisson_process.mat';...
    'gene_expression','toggle_switch.m';...
    'gene_expression','resonance.m';};

% Also include their initial conditions.
Init_Conds = {'[100,1,0]';...
    '[10,5,0]';...
    '[40,40,0]';...
    '[0,0,0]';...
    '[30,5,0]';...
    '[1,100,0]';};

% And their required integration time:
Print_times = {'[0:0.1:10]';...
    '[0:0.1:10]';...
    '[0:0.1:100]';...
    '[0:0.1:100]';...
    '[0:10:15000]';...
    '[0,10,2000]'};

Species_2fit={'[1 2 3]';...
    '[1 2]';...
    '[1 2]';...
    '[1]';...
    '[1 2]';...
    '[1 2]';};
%% Run SSA and save data
nmods = size(Models,1);
for i_mod = 1:nmods
 % Navigate to and load model
        A.ModelUsePresetExampleTypeDropDown.Value = Models{i_mod,1}; mod_change_d(A,[]);
        A.ModelDropDown.Value = Models{i_mod,2}; mod_change(A,[]);

        % Set SSA initial condition and print times
        A.SsaInitCondField.Value = Init_Conds{i_mod};
        A.PrintTimesEditField.Value = Print_times{i_mod};
        
        % Set number of simulation
        A.SsaNumSimField.Value = 10;
        
        % Run SSA
        run_SSA(A,[]);
        
        % Export data to Excel file
        data_save(A,[]);
        
        % Loads data
        load_data(A,[]);
        
        % Checks for species to fit
        if length(eval(Species_2fit{i_mod}))== 1
            A.ParEstX1CheckBox.Value = 1;
            A.ParEstX1DropDown.Value = 'x1'; 
            A.ParEstX2CheckBox.Value = 0;
            A.ParEstX2DropDown.Value = 'ignore';   
            A.ParEstX3CheckBox.Value = 0;
            A.ParEstX3DropDown.Value = 'ignore';
        elseif length(eval(Species_2fit{i_mod}))== 2
            A.ParEstX1CheckBox.Value = 1;
            A.ParEstX1DropDown.Value = 'x1';
            A.ParEstX2CheckBox.Value = 1;
            A.ParEstX2DropDown.Value = 'x2';    
            A.ParEstX3CheckBox.Value = 0;
            A.ParEstX3DropDown.Value = 'ignore';            
        else
            A.ParEstX1CheckBox.Value = 1;
            A.ParEstX1DropDown.Value = 'x1'; 
            A.ParEstX2CheckBox.Value = 1;
            A.ParEstX2DropDown.Value = 'x2';    
            A.ParEstX3CheckBox.Value = 1;
            A.ParEstX3DropDown.Value = 'x3';            
        end
        
        % Sets the parameter fitting method, number of iterations 
        % and runs Fit
        A.FittingAlgorithmDropDown.Value = 'fminsearch';
        A.NumIterartionsEditField.Value = 100;
        run_fit(A,[]);
        
        
end
