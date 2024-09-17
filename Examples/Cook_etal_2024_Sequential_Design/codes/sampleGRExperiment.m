function [simData,csvFile,allDataSoFar,maxAvailable] = sampleGRExperiment(dataFileName,times,Dex,NCells,iExpt,fname)

if iExpt == 1
    % Randomize data order
    Xtable = readtable(dataFileName);
    J = randperm(size(Xtable,1));
    Xtable = Xtable(J,:);
    writetable(Xtable,[dataFileName(1:end-4),'_',fname,'_permuted.csv']);
end
dataFileName = [dataFileName(1:end-4),'_',fname,'_permuted.csv'];

% load data
X = importdata(dataFileName);
Xtable = readtable(dataFileName);
indTime=find(ismember(X.textdata,'time'));
indDex=find(ismember(X.textdata,'Dex_Conc'));


% Find number of cells available for each experiment
maxAvailable = zeros(1,length(Dex)*length(times));
for j = 1:length(Dex)
    for i = 1:length(times)
        maxAvailable(i+(j-1)*length(times)) = ...
            sum(((X.data(:,indTime) == times(i)).*(X.data(:,indDex) == Dex(j))));
    end
end

% Pick 'NCells' random cells for simmulated experiment
for j = 1:length(Dex)
    simData = [];
    for i = 1:length(times)
        c = find((X.data(:,indTime) == times(i)).*...
            (X.data(:,indDex) == Dex(j)));
        if isempty(c)
            error(['No Cells left at timepoint t = ',num2str(times(i))])
        else
            if length(c) < NCells(j,i)
                error(['Not enough cells in experiment for measuring ',num2str(NCells(j,i)), ' cells at time t = ',num2str(timeMatrix(i))])
            else
                %                 ki = randperm(length(c),NCells(i)); % Picks the cells at random from the dataset
                ki = 1:NCells(j,i); % Picks the cells in order from the dataset
                c_pick = c(ki,:);
                simData = cat(1,simData,c_pick);
            end
        end
    end

    % Save simmulated experiment results
    Xnew=Xtable(simData,:);
    f = ['ExampleData/proposedExp_',fname,'_r',num2str(j),'.csv'];
    writetable(Xnew,f);
    csvFile{j} = f;

    if sum(NCells(j,:))>0
        allDataSoFar{j} = NCells(j,:);
    end

end

end