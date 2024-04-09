function [simData,csvFile] = sampleExperimentSim(dataFileName,timeMatrix,NCells)

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
%                 ki = randperm(length(c),NCells(i)); % Picks the cells at random from the dataset
                ki = 1:NCells(i); % Picks the cells in order from the dataset
                c_pick = c(ki,:);
                simData = cat(1,simData,c_pick);
            end
        end
    end
    
    % Save simmulated experiment results
    n=1;
    while exist(['purposeExp','_v',num2str(n),'.csv'])
        n=n+1;
    end
    Xnew=Xtable(simData,:);
    writetable(Xnew,['purposeExp','_v',num2str(n),'.csv']);
    csvFile = ['purposeExp','_v',num2str(n),'.csv'];

end