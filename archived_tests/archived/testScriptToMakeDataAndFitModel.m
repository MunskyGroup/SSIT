% This script fits the models to data created by the SSA tab. Metropolis
% Hastings and fminsearch are the two fit methods used. The time it takes
% for the data to be fit to and the parameters from each fit method are 
% recorded. 
%
% Before running this function, makeSsaData.m needs to be run. Once the
% data is created, this function can be run. makeSsaData.m only needs to be
% run once.
%%
close all force
clear all
 cd ..

src2app;
A = SSIT_GUI;  % Start up the app.
mod_change_d = get(eval(['A.ModelUsePresetExampleTypeDropDown']),'ValueChangedFcn');
mod_change = get(eval(['A.ModelDropDown']),'ValueChangedFcn');


%% changing models

Models = {'gene_expression','poisson_process.mat';...          
          'gene_expression','Gene_Toggle_Model.m';...
          'gene_expression','resonance.m'};%...
%           'epidemiology','sir_model_r0.mat'}; 
% Epidemiology model runs into an error with the ode function
      
% excel={'poisson';...       
%        'Gene_Toggle_Model';...
%        'resonance';...
%        'sir_model_r0'};
   
fileModel={'poisson.xlsx';...           
           'Gene_Toggle_Model.xlsx';...
           'resonance.xlsx'};%...
%            'sir_model_r0.xlsx'};
       
% Length of the vector determines how many species will be fit to
Species_2fit={'[1]';...              
              '[1 2]';...
              '[1 2]'};%...
%               '[1 2]'};
          
printTimes = {'[0:0.1:5]';...                             
               '[0:0.1:5]';...
               '[0:40:2000]'};%...
%                '[0:1:100]'};
           
Init_Conds = {'[0,0,0]';...
              '[10,5,0]';...             
              '[1,100,0]'};%...
%               '[100,1,0]'};
         
fitMethod={'fminsearch';          
           'Metropolis Hastings';};
 %% Loading and fitting different models  
fitTime=NaN*ones(length(fileModel),length(fitMethod));
update=NaN*ones(length(fileModel),length(fitMethod));
actualParm=num2cell(NaN*ones(length(fileModel)));

for i_mod=1%:length(Models(:,1))
       A.ModelUsePresetExampleTypeDropDown.Value = Models{i_mod,1}; mod_change_d(A,[]);
       A.ModelDropDown.Value = Models{i_mod,2}; mod_change(A,[])
       
       % Initial Conditions
%        A.SsaInitCondField.Value = Init_Conds{i_mod};
       
        % Set number of simulation
%         A.SsaNumSimField.Value = 100;
        
       % run SSA
%        A.PrintTimesEditField.Value = printTimes{i_mod};
%        A.SsaRunButton.ButtonPushedFcn(A, [])
       % you can ignore the error that pops up here for now.
       
       % save results
%        exportSsaData(A,fileModel{i_mod});
       %% open SSA results in data fitting tab
       fn = ['tests/test_data/' fileModel{i_mod}];
       pn = [];
       loadDataFromFile(A,fn,pn);
       
       % Choses the parameters to fit to
        if length(eval(Species_2fit{i_mod}))== 1
            A.ParEstX1CheckBox.Value = 1;
            A.ParEstX1DropDown.Value = {'x1'}; 
            A.ParEstX2CheckBox.Value = 0;
            A.ParEstX2DropDown.Value = {'ignore'};   
            A.ParEstX3CheckBox.Value = 0;
            A.ParEstX3DropDown.Value = {'ignore'};
        elseif length(eval(Species_2fit{i_mod}))== 2
            A.ParEstX1CheckBox.Value = 1;
            A.ParEstX1DropDown.Value = {'x1'};
            A.ParEstX2CheckBox.Value = 1;
            A.ParEstX2DropDown.Value = {'x2'};    
            A.ParEstX3CheckBox.Value = 0;
            A.ParEstX3DropDown.Value = {'ignore'};            
        else
            A.ParEstX1CheckBox.Value = 1;
            A.ParEstX1DropDown.Value = {'x1'}; 
            A.ParEstX2CheckBox.Value = 1;
            A.ParEstX2DropDown.Value = {'x2'};    
            A.ParEstX3CheckBox.Value = 1;
            A.ParEstX3DropDown.Value = {'x3'};            
        end
       
       A.ParEstFitTimesList.Value = A.ParEstFitTimesList.Items;
       userSelectCondition(A)
 
       actualParm = A.fit_parameters_table.Data 
       % Parameter values used to generate the Stochastic data
     
       %% change to ham. mc.
       A.FspErrorTolField.Value = 0.05;
      
              
       for i_fit=1:length(fitMethod)
           A.FittingAlgorithmDropDown.Value = fitMethod{i_fit};

           tFit= tic;
           ssit.parest.fitModel2Data(A);
           ssit.parest.updateModelSolveAndCompareToData(A);
           fitTime(i_mod,i_fit)=toc(tFit);
           % time for 'fit model button to run'

%            tic
%            ssit.parest.updateModelSolveAndCompareToData(A);
%            update(i_mod,i_fit)=toc 
%            % time for 'Update Model and Solve' button to run
           
           parmfitMethod{i_fit} = A.fit_parameters_table.Data;   
       end
       
       parameterValues{i_mod} = {actualParm;
                          parmfitMethod{1};% fiminsearch Parameters
                          parmfitMethod{2}}% Metropolis Hastings Parameters

end
