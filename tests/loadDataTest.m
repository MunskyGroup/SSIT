%% Test the loadData function by setting conditions on Eric Ron's GR data 
%% to filter by 3 concentrations of Dex: 1,10,100nM
close all 
clear
addpath(genpath('../src'));

    % Create SSIT model.
    ModelGR = SSIT;
    ModelGR.species = {'cytGR';'nucGR'};
    ModelGR.initialCondition = [20;1];
    ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';...
        'kg1';'gg1*cytGR';'gg2*nucGR'};
    ModelGR.stoichiometry = [-1,1,1,-1,0;...
        1,-1,0,0,-1];
    ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
        'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',14e-5;...
        'gg1',1e-5;'gg2',1e-6;'MDex',5;'Dex0',100});
    ModelGR.summarizeModel    
    ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};

    % Associate GR data with different instances of model (1,10,100nm Dex)
    GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
        '10','10',102,'GR Fit (10nM Dex)';...
        '100','100',103,'GR Fit (100nM Dex)'};
    ModelGRfit = cell(1,size(GRfitCases,1));
    
    % Loop to load data according to Dex concentration
    for i=1:3
        ModelGRfit{i} = ModelGR.loadData("test_data/GR.csv",...
                                        {'nucGR','normGRnuc';'cytGR','normGRcyt'},...
                                        {'dex_conc',GRfitCases{i,1}});  
        % Print a few columns of the data, including Dex conc. to see if
        % the filtering (loadData conditions worked)
        disp(ModelGRfit{i}.dataSet.DATA(1:3,4:6))
    end
