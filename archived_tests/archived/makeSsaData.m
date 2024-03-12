% This function creates the stochastic data that will be fit to. This
% function needs to be run before tetsScriptToMakeDataAndFitModel.m
%%
close all force
clear all
% cd ..

src2app;
A = SSIT_GUI;  % Start up the app.
mod_change_d = get(eval(['A.ModelUsePresetExampleTypeDropDown']),'ValueChangedFcn');
mod_change = get(eval(['A.ModelDropDown']),'ValueChangedFcn');


%% changing models

Models = {'gene_expression','poisson_process.mat';...          
          'gene_expression','Gene_Toggle_Model.m';...
          'gene_expression','resonance.m'};      
   
fileModel={'poisson.xlsx';...           
           'Gene_Toggle_Model.xlsx';...
           'resonance.xlsx'};
          
printTimes = {'[0:0.1:5]';...                             
               '[0:0.1:5]';...
               '[0:40:2000]'};
           
Init_Conds = {'[0,0,0]';...
              '[10,5,0]';...             
              '[1,100,0]'};
          
for i_mod=1:length(Models(:,1))
       A.ModelUsePresetExampleTypeDropDown.Value = Models{i_mod,1}; mod_change_d(A,[]);
       A.ModelDropDown.Value = Models{i_mod,2}; mod_change(A,[])
       
       % Initial Conditions
       A.SsaInitCondField.Value = Init_Conds{i_mod};
       
        % Set number of simulation
        A.SsaNumSimField.Value = 100;
        
       % run SSA
       A.PrintTimesEditField.Value = printTimes{i_mod};
       A.SsaRunButton.ButtonPushedFcn(A, [])
       % you can ignore the error that pops up here for now.
       
       % save results
       exportSsaData(A,fileModel{i_mod});
end