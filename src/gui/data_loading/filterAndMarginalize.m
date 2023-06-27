function app = filterAndMarginalize(~,~,app)
%This function applies the marginalization and filtering to the data
%structure defined by the user and stores the edited data to be plotted.

% colHeaders = app.DataLoadingAndFittingTabOutputs.columns;

%% Find relevant Indicies
histDataRaw = app.DataLoadingAndFittingTabOutputs.dataTable;
[rowNum,~] = size(histDataRaw);
% [marginalMatRow,~] = size(app.DataLoadingAndFittingTabOutputs.marginalMatrix);

Nd = length(app.SpeciesForFitPlot.Items);

% Indices for columns that should be ignored in data.
timeInd =  find(app.DataLoadingAndFittingTabOutputs.marginalMatrix(Nd+1,:));

% Record numbers of observations (last column) for each state (first N-1
% columns).  If no counts are given, assume one for each observation.
% Index corresponding to counts or probabilities for specific state.
countsInd = find(app.DataLoadingAndFittingTabOutputs.marginalMatrix(Nd+2,:));
if isempty(countsInd)                          % if no counts selected by user adds an arrary of 1s to last ind
    histDataRaw(:,end+1) = {1};
    valsHistDataInd = size(histDataRaw,2);
else
    valsHistDataInd = countsInd;                                        % find ind for values in spTensor when user selects counts from popout
end

%% Filter on Specified Conditions
% This should remove rows of the data that do not correspond to the logical
% conditions required by the user.  
% This seems to be working for equality conditions (e.g., nRep = '0'), but 
% not yet for more complicated conditions -- will require future testing.

% this part determines which variables are conditioned over
histDataStr = cellfun(@num2str,histDataRaw,'un',0);                     % convert double in cell to str
if ~isa(app.DataLoadingAndFittingTabOutputs.conditionOnArray,'double')                                      % Check if user wants to filter for specified conditions
    loopCond = cellfun(@(s) ~contains('0',s),app.DataLoadingAndFittingTabOutputs.conditionOnArray);
    loopCond = find(loopCond);
else
    loopCond = 0;
end

% this part finds the data that obeys all conditions and removes the rest
% from the data set.
if loopCond
    for iC = 1:length(loopCond)
        condStr = app.DataLoadingAndFittingTabOutputs.conditionOnArray(loopCond(iC));                       % find string to condition on
        condCell = cellstr(condStr);                                    % convert string to cell
        condIndLoc = ismember(histDataStr(:,loopCond(iC)),condCell)';
        histDataStr = histDataStr(condIndLoc,:);                        % only keep idx rows
        clearvars condInd rowNum
        [rowNum,~] = size(histDataStr);
    end
end

%% Reduce loaded and filtered data to remove marginalized quantities.
% Find indices which not marginalized over and which do not contain count
% information.

% Find indices of data corresponding to x1, x2, x3, ...
app.DataLoadingAndFittingTabOutputs.xInd = zeros(1,Nd);
for i=1:Nd
    J = find(app.DataLoadingAndFittingTabOutputs.marginalMatrix(i,:));
    if ~isempty(J)
        app.DataLoadingAndFittingTabOutputs.xInd(i) = J;
    end
end
x1x2x3Ind = app.DataLoadingAndFittingTabOutputs.xInd;
x1x2x3Ind = x1x2x3Ind(x1x2x3Ind~=0);                                    % get rid of zero indicies if there are < 3 species

histDataStr = histDataStr(:,[timeInd,x1x2x3Ind,valsHistDataInd]);

%% Map values in data structure to unique key indentifiers
% This section makes a list of all unique states that are included in the
% data.
%Create Map
loopVec = [1:size(histDataStr,2)];
margMap = {};  % initialize
histDataStr = cellfun(@str2num,histDataStr,'un',0); % convert data back to numerical values to manipulate
% histDataStrMAT = cell2mat(histDataStr);             % convert to matrix

for i = 1:length(loopVec)
    vecData = [histDataStr{:,loopVec(i)}];
    uniqueVals = unique(vecData);
    if i==1||~isnumeric(uniqueVals)||~min(floor(uniqueVals)==uniqueVals)||min(uniqueVals)<0
        margKeySet = [0:length(uniqueVals)-1];
    else
        margKeySet = uniqueVals;
    end
    if length(uniqueVals)==1
        uniqueVals = uniqueVals(1);
    end
    margMap = [margMap, {containers.Map(uniqueVals,margKeySet)}];
end
%Convert Mapped properties in data structure
for iR = 1:rowNum
    for iL = 1:length(loopVec)
        histDataStr(iR,loopVec(iL)) = values(margMap{iL}, histDataStr(iR,loopVec(iL)));
    end
end
%% Create sparse tensor
% histDataStr = cellfun(@num2str,histDataStr,'un',0);                     % convert data to strings
% histDataStr = cellfun(@str2num,histDataStr,'un',0);                     % convert data back to numerical values to manipulate
histDataStr = cell2mat(histDataStr);                                    % convert to matrix to use sptensor
subsHistData = histDataStr(:,1:end-1);
valsHistData = histDataStr(:,end);
subsHistData = subsHistData+1;            % bump up species counts so ind starts from 1
app.DataLoadingAndFittingTabOutputs.dataTensor = sptensor(subsHistData, valsHistData);                     % Create Sparse Tensor
% lrgInd = find(app.DataLoadingAndFittingTabOutputs.collapseInd>valsHistDataInd);
% app.DataLoadingAndFittingTabOutputs.collapseInd(lrgInd) = app.DataLoadingAndFittingTabOutputs.collapseInd(lrgInd) - 1;
close(gcf)
end