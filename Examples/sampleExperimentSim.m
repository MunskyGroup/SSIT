function [simData,csvFile] = sampleExperimentSim(dataFileName,timeMatrix,NCells)
% sampleExperimentSim -- Generates csv files from a real data set to create
%                        a smaller data subset and saves it as a csv file.
% arguments:
%
%   dataFileName -- The location of the experimental data set given as a
%                   string.
%                   Ex: dataFileName = '../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv'
%
%   timeMatrix --   A 1 by N vector of the experiment's measurement times.
%                   Ex: timeMatrix = [0 10 20 30 40 50 60 75 90 120 150 180]
%
%   NCells --       A vector the same size as the input 'timeMatrix' that
%                   indicates how many cells will be measured at each time
%                   point. The nth value of 'NCells' is the number of cells
%                   measured at the nth time point.
%                   Ex. In the above example for timeMatrix if
%                   NCells = [100 0 0 0 0 0 50 0 0 0 0 200], 100 cella at
%                   t = 0 min, 50 cells at t = 60 min and 200 cells at 180
%                   min are picked
%
% Outputs:
%    
%   csvFile -- A subsample of the selected cells are saved into a csv file
%              with the name saved as the variable "csvFile".
%
%   simData -- A vector of the selected rows from the original data set.
%              This shows which cells were selected from the original data
%              set.


    % load data
    X = importdata(dataFileName);
    Xtable = readtable(dataFileName);
    ind=find(ismember(X.textdata,'time_index'));
    
    % Pick 'NCells' random cells for simmulated experiment
    simData = [];
    for i = 1:length(timeMatrix) 
        c = find(X.data(:,ind) == timeMatrix(i));  
        if isempty(c)
            error(['No Cells at timepoint t = ',num2str(timeMatrix(i))])
        else
            if length(c) < NCells(i)
                error(['Not enough cells in experiment for measuring ',num2str(NCells(i)), ' cells at time t = ',num2str(timeMatrix(i))])
            else
                ki = randperm(length(c),NCells(i)); % Picks the cells at random from the dataset
                % ki = 1:NCells(i); % Picks the cells in order from the dataset
                c_pick = c(ki,:);
                simData = cat(1,simData,c_pick);
            end
        end
    end
    
    % Save simmulated experiment results
    n=1;
    while exist(['purposeExp','_v',num2str(n),'.csv'],'file')
        n=n+1;
    end
    Xnew=Xtable(simData,:);
    writetable(Xnew,['purposeExp','_v',num2str(n),'.csv']);
    csvFile = ['purposeExp','_v',num2str(n),'.csv'];

end