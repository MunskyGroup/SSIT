function testResults = install(runTests,runExamples,overwritePropFuns,saveSearchPath,publishHtml,isCluster)
arguments
    runTests = false
    runExamples = false
    overwritePropFuns = []
    saveSearchPath = []
    publishHtml = false
    isCluster = false
end

if isCluster
    % runTests = false;
    % Turn off figures.
    setenv('QT_QPA_PLATFORM','offscreen');
    set(0,'DefaultFigureVisible','off')
    runExamples = false;
    tests2run = {'poissonTest','poisson2Dtest','poissonTVtest',...
        'miscelaneousTests','multiModelTests','modelReductionTest'};
else
    tests2run = {'poissonTest','poisson2Dtest','poissonTVtest',...
        'miscelaneousTests','multiModelTests','modelReductionTest','testGui'};
end
    
% Check that 
weAreIn = pwd;
J = strfind(weAreIn,filesep);
if ~strcmpi(weAreIn(J(end)+1:end),'SSIT')
    if isCluster
        error('Not in correct SSIT directory -- cannot install')
    else
        abort = questdlg({'You appear not to be in the SSIT Directory (or you changed its name).';'Your current directory is';weAreIn;'Do you wish to abort?'}, ...
            'Confirm Action', ...
            'Yes','No','Yes');
        if strcmpi(abort,'Yes')
            return
        end
    end
end

% Set the path to include all SSIT codes.  
addpath(weAreIn)
addpath(genpath('src'));

% Install MEX codes
disp('Installing MEX codes for Expokit')
try
    ssit.fsp_ode_solvers.build_expokit;
    disp('MEX code for Expokit installed')   
catch ME
    disp('MEX code installation failure. Analyses will default to native MATLAB version.')
    disp(['Error details: ' ME.message])
    if ~isempty(ME.cause)
        for k = 1:numel(ME.cause)
            disp(['Cause: ' ME.cause{k}.message])
        end
    end
end

if ~exist("tmpPropensityFunctions","dir")
    disp('Creating director "tmpPropensityFunctions".')
    mkdir("tmpPropensityFunctions")
elseif ~isempty(dir("tmpPropensityFunctions"))&&~isCluster
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

% Test command line installation
try
    A1 = SSIT; clear A1
    disp('SSIT Command Tools are available.')
catch me
    me
    return
end

% Test GUI installation (if not on cluster)
if ~isCluster
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
else
    savepath 'src/pathdef.m';
end

if runTests
    % Run Tests
    origDir = pwd;              % save current directory
    cleanupTest = onCleanup(@() cd(origDir));  % guarantee return
    cd('tests')
    set(0, 'DefaultFigureVisible', 'off');
    disp('Starting Tests....')
    testResults.tests = runtests(tests2run);
    set(0, 'DefaultFigureVisible', 'on');
    clear cleanupTest
else
    testResults.tests =[];
end

if ~isCluster&&runExamples
    origDir = pwd;              % save current directory
    cleanupEx = onCleanup(@() cd(origDir));  % guarantee return
    cd('Examples')
    if ~exist("exampleLogs","dir")
        mkdir("exampleLogs")
    end
    ExampleFiles = {'example_00_AllManuscriptExamples.m'
        'example_01_CreateSSITModels.m'
        'example_02_SolveSSITModels_ODE.m'
        'example_03_SolveSSITModels_SSA.m'
        'example_04_SolveSSITModels_FSP.m'
        'example_05_SolveSSITModels_EscapeTimes.m'
        'example_06_SensitivityAnalysis.m'
        'example_07_FIM.m'
        'example_08_FIM_ExperimentDesign.m'
        'example_09_LoadingandFittingData_DataLoading.m'
        'example_10_LoadingandFittingData_MLE.m'
        'example_11_LoadingandFittingData_MH.m'
        'example_11b_LoadingandFittingData_MH_with_FIM.m'
        'example_12_ComplexModels_PDO.m'
        'example_13_ComplexModels_MultiModel.m'
        'example_14_PipelinesAndClusterComputing.m'
        'example_Benchmarking.m'
        'example_ModelReduction_benchmark.m'
        'example_SI_ABC.m'
        'example_SI_CrossValidation.m'
        'example_SI_Epidemics.m'
        'example_SI_Hybrid.m'
        'example_SI_MAPK.m'
        'example_SI_ModelReduction.m'
        'example_SI_Moments.m'
        'example_SI_MultiModel.m'
        'example_SI_SBML.m'
        'example_scRNAseq_1_BatchFitManyGenes.m'
        'example_scRNAseq_3_BatchFitManyGenes.m'         };

    completed = zeros(1,length(ExampleFiles),'logical');
    for iEx = 1:length(ExampleFiles)
        try 
            if publishHtml
                set(0, 'DefaultFigureVisible', 'on');
                publish(ExampleFiles{iEx})
                close all
            else
                set(0, 'DefaultFigureVisible', 'off');
                tic; out = evalc(ExampleFiles{iEx}); timeToc = toc;
                fid = fopen(['exampleLogs/output',ExampleFiles{iEx},'.txt'],'w');
                fprintf(fid,'%s', out);
                fclose(fid);
                completed(iEx) = true;
                disp([ExampleFiles{iEx},' succeeded in ',num2str(timeToc),' s; logfile in ','exampleLogs/output',ExampleFiles{iEx},'.txt'])
            end
            close all
        catch me
            set(0, 'DefaultFigureVisible', 'on');
            disp([ExampleFiles{iEx},' failed with message:'])
            me
            close all
        end
    end
    clear cleanupEx
else
    testResults.examplesCompleted = [];
end
