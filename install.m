function testResults = install(runTests,runExamples)
arguments
    runTests = false
    runExamples = false
end

% Set the path to include all SSIT codes.  
addpath(genpath('src'));

try
    A1 = SSIT; clear A1
    disp('SSIT Command Tools are available.')
catch me
    me
end

try
    A2 = SSITGUI; close(A2.UIFigure);
    disp('SSIT Command Tools and SSIT GUI are available.')
catch me
    me
end

% Run Tests
if runTests
    origDir = pwd;              % save current directory
    cleanupTest = onCleanup(@() cd(origDir));  % guarantee return
    cd('tests')
    testResults.tests = runtests({'poissonTest','poisson2Dtest','poissonTVtest',...
        'miscelaneousTests','multiModelTests','modelReductionTest','testGui'});
    clear cleanupTest
else
    testResults.tests =[];
end


if runExamples
    origDir = pwd;              % save current directory
    cleanupEx = onCleanup(@() cd(origDir));  % guarantee return
    cd('Examples')
    ExampleFiles = {'example_1_CreateSSITModels'
        'example_2_SolveSSITModels_ODE'
        'example_3_SolveSSITModels_SSA'
        'example_4_SolveSSITModels_FSP'
        'example_5_SolveSSITModels_EscapeTimes'
        'example_6_SensitivityAnalysis'
        'example_7_FIM'
        'example_8_LoadingandFittingData_DataLoading'
        'example_8b_LoadingandFittingData_SimulatingData'
        'example_9_LoadingandFittingData_MLE'
        'example_10_LoadingandFittingData_MH'
        'example_10b_LoadingandFittingData_MH_with_FIM'
        'example_11_LoadingandFittingData_CrossValidation'
        'example_12_ComplexModels_MultiModel'
        'example_13_ComplexModels_Hybrid'
        'example_14_ComplexModels_ModelReduction'
        'example_15_PDO'
        'example_16_PipelinesAndClusterComputing'
        'example_17_ABC'
        'example_18_Moments'};
    completed = zeros(1,length(ExampleFiles),'logical');
    for iEx = 1:length(ExampleFiles)
        try         
            tic; out = evalc(ExampleFiles{iEx}); timeToc = toc;
            fid = fopen(['exampleLogs/output',ExampleFiles{iEx},'.txt'],'w');
            fprintf(fid,'%s', out);
            fclose(fid);
            completed(iEx) = true;
            disp([ExampleFiles{iEx},' succeeded in ',num2str(timeToc),' s; logfile in ','exampleLogs/output',ExampleFiles{iEx},'.txt'])
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
