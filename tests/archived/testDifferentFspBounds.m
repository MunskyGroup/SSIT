%% Start app and create links to call back functions
close all force
clear all
cd ..
A = SSIT_GUI;  % Start up the app.

mod_change_d = get(eval(['A.ModelUsePresetExampleTypeDropDown']),'ValueChangedFcn');
mod_change = get(eval(['A.ModelDropDown']),'ValueChangedFcn');
run_fsp = get(eval('A.FspRunButtom'),'ButtonPushedFcn'); 

%% 1) Make list of models to test
Models = {'epidemiology','sir_model_r0.mat';...
    'gene_expression','Gene_Toggle_Model.m';...
    'social','opinion_formation.mat';...
    'gene_expression','poisson_process.m';...
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

%% 2) Make list of FSP bounds to test
FspTabOutputs.bounds(1) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3'}};
FspTabOutputs.bounds(2) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1-3)*(x2-3)'}};
FspTabOutputs.bounds(3) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'x1+x2+x3'}};
FspTabOutputs.bounds(4) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1-3)*(x2-3)';'x1+x2+x3'}};
FspTabOutputs.bounds(5) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(2*x1-3)^2*(x2-3)'}};
FspTabOutputs.bounds(6) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'x1^2+x2^2+x3^2'}};
FspTabOutputs.bounds(7) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1^2+x2^2+x3^2)^(1/2)'}};
FspTabOutputs.bounds(8) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'x1*x2/2'}};
FspTabOutputs.bounds(9) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'x3*x2/2'}};
FspTabOutputs.bounds(10) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'x1*x3/2'}};
FspTabOutputs.bounds(11) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1*x2/2)^2'}};
FspTabOutputs.bounds(12) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x3*x2/2)^2'}};
FspTabOutputs.bounds(13) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1*x3/2)^2'}};
FspTabOutputs.bounds(14) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'-x1^2-x2-x3'}};
FspTabOutputs.bounds(15) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'-x1-x2^2-x3'}};
FspTabOutputs.bounds(16) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'-x1-x2-x3^2'}};
FspTabOutputs.bounds(17) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'x1+x2+x3';'(x1+x2+x3)*2'}};
FspTabOutputs.bounds(18) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'-(x1+x2+x3)';'x1'}};
FspTabOutputs.bounds(19) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'-(x1+x2+x3)';'x2'}};
FspTabOutputs.bounds(20) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'-(x1+x2+x3)';'x2'}};
FspTabOutputs.bounds(21) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1+x2+x3)';'x1'}};
FspTabOutputs.bounds(22) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1+x2+x3)';'x2'}};
FspTabOutputs.bounds(23) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1+x2+x3)';'x2'}};
FspTabOutputs.bounds(24) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1+x2+x3)';'x1*x2'}};
FspTabOutputs.bounds(25) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1+x2+x3)';'x3*x2'}};
FspTabOutputs.bounds(26) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1+x2+x3)';'x1*x3'}};
FspTabOutputs.bounds(27) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1+x2+x3)';'x1*x2*x3'}};
FspTabOutputs.bounds(28) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1/50)^2+(x2/50)^2';'x1*x2*x3'}};
FspTabOutputs.bounds(29) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x1/50)^2+(x3/50)^2';'x1*x2*x3'}};
FspTabOutputs.bounds(30) = {{'-x1';'-x2';'-x3';'x1';'x2';'x3';'(x3/50)^2+(x2/50)^2';'x1*x2*x3'}};

%% 3) Iterate to test different bounds on different models
nbnds = length(FspTabOutputs.bounds);
nmods = size(Models,1);
Compute_Time1 = NaN*ones(nmods+1,nbnds+1);
Compute_Time_Fin1 = NaN*ones(nmods+1,nbnds+1);
FSP_Size1 = NaN*ones(nmods+1,nbnds+1);

Compute_Time=num2cell(Compute_Time1);
Compute_Time_Fin=num2cell(Compute_Time_Fin1);
FSP_Size=num2cell(FSP_Size1);

for i_mod = 1:nmods
        
        Compute_Time{i_mod+1,1}=Models{i_mod,2}
        Compute_Time_Fin{i_mod+1,1}=Models{i_mod,2}
        FSP_Size{i_mod+1,1}=Models{i_mod,2}
        
        for i_bnds = 1:nbnds
            
            FSP_Size{1,i_bnds+1}=Bounds{i_bnds}{end}
            Compute_Time_Fin{1,i_bnds+1}=Bounds{i_bnds}{end}
            Compute_Time{1,i_bnds+1}=Bounds{i_bnds}{end}
            
        % Navigate to and load model
        A.ModelUsePresetExampleTypeDropDown.Value = Models{i_mod,1}; mod_change_d(A,[]);
        A.ModelDropDown.Value = Models{i_mod,2}; mod_change(A,[]);

        % Set FSP iniitial condition and print times
        A.FspInitCondField.Value = Init_Conds{i_mod};
        A.FspPrintTimesField.Value = Print_times{i_mod};
       
        % Set FSP bound functions
        A.FspConstraintTable.Data = {};
        A.FspConstraintTable.Data(1:length(FspTabOutputs.bounds{i_bnds}),1) = FspTabOutputs.bounds{i_bnds};
        A.FspConstraintTable.Data(1:length(FspTabOutputs.bounds{i_bnds}),2) = {'<'};
        A.FspConstraintTable.Data(1:3,3) = {0};        
        A.FspConstraintTable.Data(4:length(FspTabOutputs.bounds{i_bnds}),3) = {1};        
        
        % compute time to expand and find FSP solution
        tic
<<<<<<< HEAD
        run_fsp(A,[]);
        Compute_Time(i_mod,i_bnds) = toc
        FSP_Size(i_mod,i_bnds) = length(A.FspTabOutputs.solutions{end}.states)
=======
        %run_fsp(A,[]);
        runFsp(A);
        Compute_Time{i_mod+1,i_bnds+1} = toc
        FSP_Size{i_mod+1,i_bnds+1} = length(A.FSP_Results_Array{end}.states)
>>>>>>> origin/main1_testbounds
        
        % compute time just for FSP solution
        tic
        run_fsp(A,[]);
        Compute_Time_Fin{i_mod+1,i_bnds+1} = toc
        
    end
end
