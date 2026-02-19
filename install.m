function testResults = install(runTests,runExamples,overwritePropFuns,saveSearchPath,publishHtml)
arguments
    runTests = false
    runExamples = false
    overwritePropFuns = []
    saveSearchPath = []
    publishHtml = false
end

% Set the path to include all SSIT codes.  
addpath(genpath('src'));

if ~exist("tmpPropensityFunctions")
    disp('Creating director "tmpPropensityFunctions".')
    mkdir("tmpPropensityFunctions")
elseif ~isempty(dir("tmpPropensityFunctions"))
    if isempty(overwritePropFuns)
        overwritePropFuns = questdlg({'Directory "tmpPropensityFunctions" already exists.','Do you wish to delete for a clean installation?'}, ...
            'Confirm Action', ...
            'Yes','No','No');
    end
    switch overwritePropFuns
        case 'Yes'
            rmdir("tmpPropensityFunctions",'s');
            mkdir("tmpPropensityFunctions");
    end
end

pathToPropensityFuns = append(pwd,filesep,'tmpPropensityFunctions');
configFile = append('src',filesep,'SSITconfig.mat');
if ~exist(configFile,'file')
    save(configFile,'pathToPropensityFuns');
else
    save(configFile,'pathToPropensityFuns','-append');
end

try
    A1 = SSIT; clear A1
    disp('SSIT Command Tools are available.')
catch me
    me
    return
end

try
    A2 = SSITGUI; close(A2.UIFigure);
    disp('SSIT Command Tools and SSIT GUI are available.')
catch me
    me
    return
end

if isempty(saveSearchPath)
    saveSearchPath = questdlg({'Installation is successful.','Do you wish to save the MATLAB search path','to use SSIT in future sessions?'}, ...
        'Confirm Action', ...
        'Yes','No','No');
end
switch saveSearchPath
    case 'Yes'
        savepath
end

if runTests
    % Run Tests
    origDir = pwd;              % save current directory
    cleanupTest = onCleanup(@() cd(origDir));  % guarantee return
    cd('tests')
    testResults.tests = runtests({'poissonTest','poisson2Dtest','poissonTVtest',...
        'miscelaneousTests','multiModelTests','modelReductionTest','testGui'})
    clear cleanupTest
else
    testResults.tests =[];
end

if runExamples
    origDir = pwd;              % save current directory
    cleanupEx = onCleanup(@() cd(origDir));  % guarantee return
    cd('Examples')
    if ~exist("exampleLogs","dir")
        mkdir("exampleLogs")
    end
    ExampleFiles = {'example_1_CreateSSITModels'
        'example_2_SolveSSITModels_ODE'
        'example_3_SolveSSITModels_SSA'
        'example_4_SolveSSITModels_FSP'
        'example_5_SolveSSITModels_EscapeTimes'
        'example_6_SensitivityAnalysis'
        'example_7_FIM'
        % 'example_8_LoadingandFittingData_DataLoading'
        % 'example_8b_LoadingandFittingData_SimulatingData'
        % 'example_9_LoadingandFittingData_MLE'
        % 'example_10_LoadingandFittingData_MH'
        % 'example_10b_LoadingandFittingData_MH_with_FIM'
        % 'example_11_LoadingandFittingData_CrossValidation'
        % 'example_12_ComplexModels_ModelReduction'
        % 'example_13_ComplexModels_Hybrid'
        % 'example_14_ComplexModels_PDO'
        % 'example_15_ComplexModels_MultiModel'
        % 'example_16_PipelinesAndClusterComputing'
        % 'example_SI_ABC'
        % 'example_SI_Moments'
        };
    completed = zeros(1,length(ExampleFiles),'logical');
    for iEx = 1:length(ExampleFiles)
        try         
            if publishHtml
                publish(ExampleFiles{iEx})
            else
                tic; out = evalc(ExampleFiles{iEx}); timeToc = toc;
                fid = fopen(['exampleLogs/output',ExampleFiles{iEx},'.txt'],'w');
                fprintf(fid,'%s', out);
                fclose(fid);
                completed(iEx) = true;
                disp([ExampleFiles{iEx},' succeeded in ',num2str(timeToc),' s; logfile in ','exampleLogs/output',ExampleFiles{iEx},'.txt'])
            end
            close all
        catch me
            disp([ExampleFiles{iEx},' failed with message:'])
            me
            close all
        end
    end
    clear cleanupEx
else
    testResults.examplesCompleted = [];
end
