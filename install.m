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
    runTests = false;
    runExamples = false;
end
    
% Check that 
weAreIn = pwd;
J = strfind(weAreIn,filesep);
if ~strcmpi(weAreIn(J(end)+1:end),'SSIT')
    if isCluster
        error('Not in correct SSIT directory -- cannot install')
    else
        abort = questdlg({'You appear not to be in the SSIT Directory.';'Your current directory is';weAreIn;'Do you wish to abort?'}, ...
            'Confirm Action', ...
            'Yes','No','Yes');
        if strcmpi(abort,'Yes')
            return
        end
    end
end

% Set the path to include all SSIT codes.  
addpath(genpath('src'));

% Install MEX codes
disp('Installing MEX codes for Expokit')
try
    ssit.fsp_ode_solvers.build_expokit;
    disp('MEX code for Expokit installed')   
catch ME
    disp('MEX code installation failure. Analyses will default to native MATLAB version.')
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
    testResults.tests = runtests({'poissonTest','poisson2Dtest','poissonTVtest',...
        'miscelaneousTests','multiModelTests','modelReductionTest','testGui'})
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
        % 'example_11_ComplexModels_PDO'
        % 'example_12_PipelinesAndClusterComputing'
        % 'example_SI_ABC'
        % 'example_SI_Moments'
        % 'example_SI_CrossValidation'
        % 'example_SI_ModelReduction'
        % 'example_SI_Hybrid'
        % 'example_SI_MultiModel'
        };
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
