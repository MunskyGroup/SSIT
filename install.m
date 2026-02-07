function testResults = install(runTests)
arguments
    runTests = false
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
    cleanupObj = onCleanup(@() cd(origDir));  % guarantee return
    cd('tests')
    testResults = runtests({'poissonTest','poisson2Dtest','poissonTVtest',...
        'miscelaneousTests','multiModelTests','modelReductionTest','testGui'});
else
    testResults =[];
end
